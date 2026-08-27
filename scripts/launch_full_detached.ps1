# =============================================================================
# scripts/launch_full_detached.ps1
#
# Launch the full-mode pipeline so it SURVIVES the terminal / IDE / agent that
# started it.
#
# WHY THIS EXISTS. Start-Process leaves the new process inside the caller's
# process tree and job object. On 2026-08-23 the Claude desktop app logged two
# MoAppHang Windows Error Reporting events at 21:18:25 and 21:18:30; the
# pipeline's worker future was interrupted at 21:18:07 and the run died with
#   Future ('diagnostics_calibrated_malawi_women_zinc') ... interrupted (pid ...)
# No R error, no allocation failure -- the workers were simply killed when the
# host app went down. A multi-hour run must not be hostage to that.
#
# Win32_Process::Create parents the new process to the WMI provider host
# instead, so it is unaffected by the launching app's lifetime. It needs no
# elevation and cannot redirect output itself, hence the cmd /c wrapper.
#
# Usage:  powershell -ExecutionPolicy Bypass -File scripts/launch_full_detached.ps1
# =============================================================================

param(
  [string]$RScript = "C:\Program Files\R\R-4.4.2\bin\x64\Rscript.exe",
  [string]$Script  = "scripts/run_full_mode.R",
  [string]$OutLog  = "pipeline_full.log",
  [string]$ErrLog  = "pipeline_full.err",
  [int]   $Workers = 0        # 0 = let run_full_mode.R size it from RAM
)

$proj = (Resolve-Path "$PSScriptRoot\..").Path
Set-Location $proj

# ---------------------------------------------------------------------------
# Stop only THIS project's R processes.
#
# A blanket `Get-Process Rscript,Rterm | Stop-Process` kills every R session on
# the machine -- another project's job, an RStudio console, someone else's
# pipeline. Match on the command line instead, which for our processes always
# references this directory, run_full_mode.R, or the callr script that the
# targets master runs.
# ---------------------------------------------------------------------------
$escaped = [regex]::Escape($proj)
$allR = Get-CimInstance Win32_Process -Filter "Name='Rscript.exe' OR Name='Rterm.exe'"
# Seeds are the processes we can positively identify from their command line.
# future/multisession WORKERS cannot be matched that way -- they run
# "parallel:::.workRSOCK()" with no reference to the project -- so they are
# picked up by walking down from a seed. Without this the cleanup silently
# missed every worker and left 4-6 orphans running.
$seedR = @($allR | Where-Object { $_.CommandLine -match $escaped -or
                                  $_.CommandLine -match 'run_full_mode|callr-scr' })
$pids = [System.Collections.Generic.HashSet[int]]::new()
$seedR | ForEach-Object { [void]$pids.Add([int]$_.ProcessId) }
for ($i = 0; $i -lt 4; $i++) {
  foreach ($q in $allR) { if ($pids.Contains([int]$q.ParentProcessId)) { [void]$pids.Add([int]$q.ProcessId) } }
}
$mine = @($allR | Where-Object { $pids.Contains([int]$_.ProcessId) })
if ($mine) {
  Write-Host "stopping $($mine.Count) stale R process(es) from this project"
  # Kill the master(s) first so they cannot re-dispatch to workers mid-teardown.
  $mine | Sort-Object { if ($_.CommandLine -match 'callr-scr') { 0 } else { 1 } } |
    ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
  Start-Sleep -Seconds 4
} else {
  Write-Host "no stale R processes from this project"
}

# Keep the previous log rather than silently overwriting the evidence.
if (Test-Path $OutLog) {
  $stamp = (Get-Item $OutLog).LastWriteTime.ToString('yyyyMMdd-HHmmss')
  Move-Item $OutLog "pipeline_full_$stamp.log" -Force
  Write-Host "previous log kept as pipeline_full_$stamp.log"
}

# MN_TARGETS_WORKER must NOT be inherited -- _targets.R uses it to tell a worker
# from the master, and a stale "1" would leave the master sequential.
Remove-Item Env:\MN_TARGETS_WORKER -ErrorAction SilentlyContinue
if ($Workers -gt 0) { $env:TARGETS_WORKERS = "$Workers" }
else { Remove-Item Env:\TARGETS_WORKERS -ErrorAction SilentlyContinue }

$cmd = 'cmd.exe /c ""{0}" "{1}" > "{2}" 2> "{3}""' -f $RScript, $Script, $OutLog, $ErrLog
$r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create `
       -Arguments @{ CommandLine = $cmd; CurrentDirectory = $proj }

if ($r.ReturnValue -ne 0) {
  Write-Error "Win32_Process::Create failed with ReturnValue=$($r.ReturnValue)"
  exit 1
}
# Hold off sleep for the duration of the run. Separate process so it can
# release the request the moment the pipeline exits, and so a failure here
# cannot take the pipeline down with it. See scripts/keep_awake.ps1.
$ka = 'powershell.exe -ExecutionPolicy Bypass -WindowStyle Hidden -File "{0}\scripts\keep_awake.ps1"' -f $proj
$kr = Invoke-CimMethod -ClassName Win32_Process -MethodName Create `
        -Arguments @{ CommandLine = $ka; CurrentDirectory = $proj }
if ($kr.ReturnValue -eq 0) {
  Write-Host "keep-awake helper started (PID $($kr.ProcessId)); releases automatically when the run ends"
} else {
  Write-Warning "keep-awake helper failed to start (ReturnValue=$($kr.ReturnValue)); the machine may sleep mid-run"
}

Write-Host "launched detached: cmd PID $($r.ProcessId) at $(Get-Date -Format 'HH:mm:ss')"
Write-Host "this process is parented to WMI, not to the calling app -- it will"
Write-Host "keep running if the terminal or agent goes away."
