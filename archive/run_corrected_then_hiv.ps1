# ===========================================================================
# run_corrected_then_hiv.ps1  (ASCII-only for Windows PowerShell 5.1)
#
# Single detached, logged, checkpointed driver:
#   1. Run the FULL corrected pipeline (tar_make, CORRECTED_SCOPE=full),
#      resume-on-restart via the _targets_corrected/ cache.
#   2. ON SUCCESS ONLY: chain HIV Stage 3 (delete its HOLD sentinel + relaunch
#      its driver). On failure/incomplete: do NOT trigger HIV; leave it held.
# Single-instance (lock file). No interactive prompts. Only ONE heavy job runs
# at a time (corrected, then HIV) - never both.
# ===========================================================================
$ErrorActionPreference = "Continue"
$root    = "C:\Users\andre\OneDrive\Documents\mn-prediction"
$rscript = "C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
$log     = Join-Path $root "corrected_full_run.log"
$lock    = Join-Path $root ".corrected_driver.lock"
$doneS   = Join-Path $root ".corrected_done"
$failS   = Join-Path $root ".corrected_failed"
$makeOut = Join-Path $root "corrected_make.out"
$makeErr = Join-Path $root "corrected_make.err"

function Log($m) {
  $ts = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
  Add-Content -Path $log -Value "$ts  $m"
}

# single instance
if (Test-Path $lock) {
  Log "Driver lock present ($lock) - another instance is running. Exiting."
  exit 0
}
Set-Content -Path $lock -Value "$PID"
Remove-Item -Path $doneS -Force -ErrorAction SilentlyContinue
Remove-Item -Path $failS -Force -ErrorAction SilentlyContinue

Log "==== CORRECTED PIPELINE DRIVER START (PID $PID) ===="
Log "Running full corrected pipeline (tar_make, CORRECTED_SCOPE=full)..."

# 1. corrected pipeline
$env:CORRECTED_SCOPE = "full"
$runScript = Join-Path $root "corrected_driver\run_corrected_full.R"
$p = Start-Process -FilePath $rscript -ArgumentList @('--vanilla', $runScript) -WorkingDirectory $root -NoNewWindow -Wait -PassThru -RedirectStandardOutput $makeOut -RedirectStandardError $makeErr
$code = $p.ExitCode
if (Test-Path $makeOut) { Get-Content $makeOut | ForEach-Object { Add-Content $log "  [make] $_" } }
if (Test-Path $makeErr) { Get-Content $makeErr | ForEach-Object { Add-Content $log "  [make:err] $_" } }
Log "tar_make process exit code: $code"

# 2. success verification
$ok = $false
if ($code -eq 0) {
  $vScript = Join-Path $root "corrected_driver\verify_corrected.R"
  $vout = Join-Path $root "corrected_verify.out"
  $v = Start-Process -FilePath $rscript -ArgumentList @('--vanilla', $vScript) -WorkingDirectory $root -NoNewWindow -Wait -PassThru -RedirectStandardOutput $vout -RedirectStandardError "$vout.err"
  if (Test-Path $vout) { Get-Content $vout | ForEach-Object { Add-Content $log "  [verify] $_" } }
  Log "verify exit code: $($v.ExitCode)"
  if ($v.ExitCode -eq 0) { $ok = $true }
} else {
  Log "tar_make exited non-zero - skipping verify."
}

# 3. completion hook
if ($ok) {
  Set-Content -Path $doneS -Value (Get-Date -Format o)
  Log "CORRECTED PIPELINE COMPLETE - SUCCESS."
  Log "==== HANDOFF: resuming HIV Stage 3 ===="
  try {
    Remove-Item -Path "C:\Users\andre\OneDrive\Documents\hiv_adherence_simulation\results\new_analyses\.stage3_done" -Force -ErrorAction SilentlyContinue
    Log "Deleted HIV .stage3_done HOLD sentinel."
    Log "Launching HIV pipeline: run_new_analyses_pipeline.ps1"
    $hivScript = "C:\Users\andre\OneDrive\Documents\hiv_adherence_simulation\run_new_analyses_pipeline.ps1"
    $hivOut = Join-Path $root "hiv_stage3_handoff.out"
    $h = Start-Process -FilePath "powershell.exe" -ArgumentList @('-NoProfile','-ExecutionPolicy','Bypass','-File', $hivScript) -NoNewWindow -Wait -PassThru -RedirectStandardOutput $hivOut -RedirectStandardError "$hivOut.err"
    Log "HIV Stage 3 driver returned (exit $($h.ExitCode)). Output -> $hivOut"
  } catch {
    Log "HIV handoff ERROR: $($_.Exception.Message)"
  }
} else {
  Set-Content -Path $failS -Value (Get-Date -Format o)
  Log "CORRECTED PIPELINE FAILED or INCOMPLETE - HIV Stage 3 NOT triggered (left held)."
}

Remove-Item -Path $lock -Force -ErrorAction SilentlyContinue
Log "==== DRIVER END ===="
