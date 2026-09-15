# ===========================================================================
# run_shap_kernel.ps1  (ASCII-only; Windows PowerShell 5.1)
#
# Detached, logged, single-instance launcher for the model-agnostic SHAP batch
# (Option B). The R driver (shap_kernel/run_shap_kernel.R) is itself
# checkpointed and resume-on-restart: re-running skips completed slices, so if
# this is killed/restarted it picks up where it left off.
#
# DOES NOT run automatically — launch it yourself per docs/SHAP_KERNEL_RUN.md.
# ===========================================================================
$ErrorActionPreference = "Continue"
$root    = "C:\Users\andre\OneDrive\Documents\mn-prediction"
$rscript = "C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
$log     = Join-Path $root "shap_kernel_run.log"
$lock    = Join-Path $root ".shap_kernel.lock"
$rout    = Join-Path $root "shap_kernel_run.out"

function Log($m) {
  $ts = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
  Add-Content -Path $log -Value "$ts  $m"
}

if (Test-Path $lock) { Log "Lock present - another SHAP run is active. Exiting."; exit 0 }
Set-Content -Path $lock -Value "$PID"

Log "==== SHAP KERNEL BATCH START (PID $PID) ===="
$script = Join-Path $root "shap_kernel\run_shap_kernel.R"
$p = Start-Process -FilePath $rscript -ArgumentList @('--vanilla', $script) `
  -WorkingDirectory $root -NoNewWindow -Wait -PassThru `
  -RedirectStandardOutput $rout -RedirectStandardError "$rout.err"
if (Test-Path $rout)        { Get-Content $rout        | ForEach-Object { Add-Content $log "  $_" } }
if (Test-Path "$rout.err")  { Get-Content "$rout.err"  | ForEach-Object { Add-Content $log "  [err] $_" } }
Log "Rscript exit code: $($p.ExitCode)"

Remove-Item -Path $lock -Force -ErrorAction SilentlyContinue
Log "==== SHAP KERNEL BATCH END ===="
