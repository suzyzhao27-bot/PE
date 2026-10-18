param(
    [Parameter(Mandatory = $true)]
    [string]$InputPath,
    [Parameter(Mandatory = $true)]
    [string]$OutputPath
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$source = (Resolve-Path -LiteralPath $InputPath).Path
$destination = [System.IO.Path]::GetFullPath($OutputPath)
if (-not (Test-Path -LiteralPath (Split-Path -Parent $destination) -PathType Container)) {
    throw 'Output parent directory must already exist.'
}

$origin = $null
try {
    $origin = New-Object -ComObject Origin.Application
    $origin.Visible = 0
    if (-not $origin.Load($source)) { throw "Origin failed to load: $source" }
    $null = $origin.Run()

    # Read saved graph X/Y source ranges. No project Save command is issued.
    $command = 'string __codex_plots$; __codex_plots$=""; doc -e LP { doc -e D { string __codex_page$=%H; string __codex_data$=%C; range -wx __codex_x = %C; range -wy __codex_y = %C; string __codex_xs$=%(__codex_x); string __codex_ys$=%(__codex_y); __codex_plots$=__codex_plots$+__codex_page$+"^"+__codex_data$+"^"+__codex_xs$+"^"+__codex_ys$+"|"; }; };'
    if (-not $origin.Execute($command)) { throw 'Origin graph-range query failed.' }
    $raw = [string]$origin.LTStr('__codex_plots$')
    $plots = @(
        foreach ($item in $raw.Split('|', [System.StringSplitOptions]::RemoveEmptyEntries)) {
            $fields = $item.Split('^')
            if ($fields.Count -ne 4) { throw "Unexpected Origin plot record: $item" }
            [ordered]@{
                graph = $fields[0]
                dataset = $fields[1]
                x_range = $fields[2]
                y_range = $fields[3]
            }
        }
    )
    if ($plots.Count -eq 0) { throw 'Origin project has no saved graph plots.' }
    $report = [ordered]@{
        source_file = [System.IO.Path]::GetFileName($source)
        method = 'Origin saved graph plot ranges; read-only COM Load, no Save'
        plots = $plots
    }
    $json = $report | ConvertTo-Json -Depth 5
    [System.IO.File]::WriteAllText($destination, $json, [System.Text.UTF8Encoding]::new($false))
    "saved_plots=$($plots.Count) output=$destination"
}
finally {
    if ($null -ne $origin) {
        try { $null = $origin.Exit() } catch { }
        [System.Runtime.InteropServices.Marshal]::FinalReleaseComObject($origin) | Out-Null
    }
}
