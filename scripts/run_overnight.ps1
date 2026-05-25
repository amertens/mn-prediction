# Robust unattended runner for the full-mode pipeline (sequential, BART kept).
# Re-launches tar_make() if it crashes; targets resumes from completed targets,
# so transient crashes (e.g. memory blips) self-heal. Stops when a run exits 0
# (all targets up to date) or after the retry cap.
$env:R_PROFILE_USER = 'NUL'
$env:PIPELINE_MODE   = 'full'
$env:TARGETS_WORKERS = '1'        # sequential; ablation future.apply uses 1 worker
$rs  = 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe'
$log = 'pipeline_full_rerun.log'
"START $(Get-Date)" | Out-File -FilePath $log -Encoding utf8
for ($i = 1; $i -le 12; $i++) {
  "=== attempt $i  $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" | Out-File -FilePath $log -Append -Encoding utf8
  if (Test-Path '_targets_full/meta/process') { Remove-Item '_targets_full/meta/process' -Force -ErrorAction SilentlyContinue }
  & $rs --vanilla -e "targets::tar_make()" *>> $log
  $code = $LASTEXITCODE
  "=== attempt $i exited with code $code ===" | Out-File -FilePath $log -Append -Encoding utf8
  if ($code -eq 0) { "ALL TARGETS COMPLETE" | Out-File -FilePath $log -Append -Encoding utf8; break }
  Start-Sleep -Seconds 30
}
"FINISHED $(Get-Date)" | Out-File -FilePath $log -Append -Encoding utf8
