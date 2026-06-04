param(
  [string]$HostName = "127.0.0.1",
  [int]$Port = 8076
)

$ErrorActionPreference = "Stop"
$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
$env:PYTHONPATH = Join-Path $repoRoot "src"
$env:TRACT_RENDERER_HOST = $HostName
$env:TRACT_RENDERER_PORT = [string]$Port

Write-Host "Starting tract-reference-renderer on ws://$HostName`:$Port"
py -3 -m tract_reference_renderer
