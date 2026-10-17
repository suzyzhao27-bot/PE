param(
    [Parameter(Mandatory = $true)]
    [string]$OutputDirectory,
    [string]$PythonExecutable
)

$ErrorActionPreference = 'Stop'
$moduleRoot = $PSScriptRoot
$repositoryRoot = (Resolve-Path -LiteralPath (Join-Path $moduleRoot '../..')).Path
if (-not $PythonExecutable) {
    $PythonExecutable = Join-Path $repositoryRoot '.venv/Scripts/python.exe'
}
if (-not (Test-Path -LiteralPath $PythonExecutable -PathType Leaf)) {
    throw "Python not found: $PythonExecutable. Specify -PythonExecutable; this script does not install dependencies."
}
$PythonExecutable = (Resolve-Path -LiteralPath $PythonExecutable).Path

$outputPath = [System.IO.Path]::GetFullPath($OutputDirectory)
$repositoryPrefix = $repositoryRoot.TrimEnd('\', '/') + [System.IO.Path]::DirectorySeparatorChar
if ($outputPath.StartsWith($repositoryPrefix, [System.StringComparison]::OrdinalIgnoreCase) -or
    $outputPath.Equals($repositoryRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw 'Output directory must be outside the repository.'
}
if (Test-Path -LiteralPath $outputPath) {
    throw "Output directory already exists; choose a fresh path: $outputPath"
}

function Invoke-ModelCommand {
    param([string]$WorkingDirectory, [string[]]$Arguments)
    $step = [ordered]@{
        working_directory = $WorkingDirectory
        command = @($PythonExecutable, '-B') + $Arguments
        exit_code = $null
    }
    Push-Location -LiteralPath $WorkingDirectory
    try {
        & $PythonExecutable -B @Arguments
        $step.exit_code = $LASTEXITCODE
        if ($LASTEXITCODE -ne 0) {
            throw "Reproduction failed: $WorkingDirectory / $($Arguments -join ' ') (exit $LASTEXITCODE)"
        }
    } finally {
        Pop-Location
        $script:runRecord.steps += $step
    }
}

$pythonVersion = (& $PythonExecutable --version 2>&1 | Out-String).Trim()
if ($LASTEXITCODE -ne 0) {
    throw "Could not query Python version: $PythonExecutable"
}
New-Item -ItemType Directory -Path $outputPath | Out-Null
$runRecord = [ordered]@{
    started_at = (Get-Date).ToString('o')
    source_directory = $moduleRoot
    output_directory = $outputPath
    python_executable = $PythonExecutable
    python_version = $pythonVersion
    status = 'in_progress'
    steps = @()
}
try {
    foreach ($part in @('core', 'fit', 'pnp')) {
        Copy-Item -LiteralPath (Join-Path $moduleRoot $part) -Destination $outputPath -Recurse
    }
    Invoke-ModelCommand (Join-Path $outputPath 'core') @('test_model.py')
    Invoke-ModelCommand (Join-Path $outputPath 'core') @('run_student.py')
    Invoke-ModelCommand (Join-Path $outputPath 'fit') @('fit_shared_overlap.py')
    Invoke-ModelCommand (Join-Path $outputPath 'pnp') @('tests_energy.py')
    Invoke-ModelCommand (Join-Path $outputPath 'pnp') @('legacy_tests.py')
    Invoke-ModelCommand (Join-Path $outputPath 'pnp') @('run_fitted_field_check.py')
    $runRecord.status = 'completed'
} catch {
    $runRecord.status = 'failed'
    throw
} finally {
    $runRecord.finished_at = (Get-Date).ToString('o')
    $runRecord | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $outputPath 'run_record.json') -Encoding UTF8
}
Write-Output "Reproduction complete; working copy and results: $outputPath"
