# =============================================================================
# scripts/keep_awake.ps1
#
# Hold a system power request for as long as this project's R processes are
# running, then release it. Nothing persists: the request lives and dies with
# this process, and no power SETTING is modified.
#
# WHY. The 2026-08-23 full run was alive and computing at 03:14 and dead by
# 05:01. The System log shows
#     05:00:08  "The system is exiting Modern Standby | Reason: Input Keyboard"
# and pipeline_full.err received a "^C" at 05:00:17. The machine slept, and the
# console process took an interrupt on wake. An interrupted LOCO target restarts
# from zero, so on a machine that sleeps nightly a multi-hour target can never
# finish no matter how many times it is relaunched.
#
# CAVEAT, stated plainly: ES_SYSTEM_REQUIRED reliably defers the classic idle
# sleep timer. On Modern Standby (S0 low-power idle) systems it is less
# absolute -- the OS can still enter DRIPS. Verify with `powercfg /requests`,
# which should list this PowerShell process under SYSTEM while it runs.
#
# Started automatically by scripts/launch_full_detached.ps1.
# =============================================================================

Add-Type -Namespace Win32 -Name Power -MemberDefinition @'
[DllImport("kernel32.dll", SetLastError = true)]
public static extern uint SetThreadExecutionState(uint esFlags);
'@

$projPattern = [regex]::Escape((Resolve-Path "$PSScriptRoot\..").Path)

$ES_CONTINUOUS       = [uint32]0x80000000
$ES_SYSTEM_REQUIRED  = [uint32]0x00000001

# ES_CONTINUOUS makes the request persist until explicitly cleared, rather than
# resetting the idle timer exactly once.
$null = [Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS -bor $ES_SYSTEM_REQUIRED)
Write-Host "keep_awake: power request held at $(Get-Date -Format 'HH:mm:ss')"

try {
  # Give the pipeline a moment to actually start before we begin checking, so we
  # do not exit during the gap between launch and the first R process.
  Start-Sleep -Seconds 60
  while ($true) {
    # Count only THIS project's R processes. Counting every Rscript on the
    # machine means an unrelated project's R session holds the machine awake
    # indefinitely after our run has finished -- which is exactly what happened
    # on 2026-08-24, when a debug_fusion.R job in another repo kept the request
    # alive after this pipeline had stopped.
    #
    # Workers carry only "parallel:::.workRSOCK()" on their command line, so
    # they are matched by ancestry from an identifiable seed rather than by text.
    $all  = Get-CimInstance Win32_Process -Filter "Name='Rscript.exe' OR Name='Rterm.exe'"
    $seed = @($all | Where-Object { $_.CommandLine -match $projPattern -or
                                    $_.CommandLine -match 'run_full_mode|callr-scr' })
    $ids = [System.Collections.Generic.HashSet[int]]::new()
    $seed | ForEach-Object { [void]$ids.Add([int]$_.ProcessId) }
    for ($i = 0; $i -lt 4; $i++) {
      foreach ($q in $all) { if ($ids.Contains([int]$q.ParentProcessId)) { [void]$ids.Add([int]$q.ProcessId) } }
    }
    if ($ids.Count -eq 0) { break }
    Start-Sleep -Seconds 30
  }
} finally {
  # Always release, including on Ctrl+C or if the loop throws. Leaving the flag
  # set would keep the machine awake indefinitely after the run is gone.
  $null = [Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS)
  Write-Host "keep_awake: power request released at $(Get-Date -Format 'HH:mm:ss')"
}
