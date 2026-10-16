param(
    [Parameter(Mandatory = $true)]
    [string[]]$InputPath,

    [Parameter(Mandatory = $true)]
    [string]$OutputRoot
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Convert-ToSafeName {
    param([string]$Name)
    $safe = $Name -replace '[^\p{L}\p{Nd}._-]+', '_'
    if ([string]::IsNullOrWhiteSpace($safe)) { return 'unnamed' }
    return $safe.Trim('_')
}

function Convert-ToCsvField {
    param($Value)
    if ($null -eq $Value) { return '""' }

    if ($Value -is [double] -or $Value -is [single] -or $Value -is [decimal]) {
        $text = $Value.ToString('R', [System.Globalization.CultureInfo]::InvariantCulture)
    }
    elseif ($Value -is [datetime]) {
        $text = $Value.ToString('o', [System.Globalization.CultureInfo]::InvariantCulture)
    }
    else {
        $text = [string]$Value
    }

    return '"' + ($text -replace '"', '""') + '"'
}

function Get-OriginSheetNames {
    param($Application)

    $command = 'string __codex_sheets$; __codex_sheets$=""; doc -e LB {string __codex_book$=%H; string __codex_sheet$=layer.name$; __codex_sheets$=__codex_sheets$+"["+__codex_book$+"]"+__codex_sheet$+"|";};'
    $null = $Application.Execute($command)
    $raw = [string]$Application.LTStr('__codex_sheets$')
    return @($raw.Split('|', [System.StringSplitOptions]::RemoveEmptyEntries))
}

function Export-OriginWorksheet {
    param(
        $Application,
        [string]$SheetRange,
        [string]$ProjectOutput
    )

    $worksheet = $Application.FindWorksheet($SheetRange)
    if ($null -eq $worksheet) {
        throw "Origin worksheet not found: $SheetRange"
    }

    $bookName = ($SheetRange -replace '^\[([^]]+)\].*$', '$1')
    $sheetName = ($SheetRange -replace '^\[[^]]+\]', '')
    $safeStem = Convert-ToSafeName "${bookName}__${sheetName}"
    $csvPath = Join-Path $ProjectOutput "$safeStem.csv"
    $metadataPath = Join-Path $ProjectOutput "$safeStem.metadata.json"

    $columns = @()
    for ($columnIndex = 0; $columnIndex -lt [int]$worksheet.Cols; $columnIndex++) {
        $column = $worksheet.Columns.Item($columnIndex)
        $columns += [ordered]@{
            index       = $columnIndex
            short_name  = [string]$column.Name
            long_name   = [string]$column.LongName
            units       = [string]$column.Units
            comments    = [string]$column.Comments
            designation = [int]$column.Type
            data_format = [int]$column.DataFormat
        }
    }

    $data = $Application.GetWorksheet($SheetRange, 0, 0, -1, -1, 0)
    $rowCount = 0
    $columnCount = 0
    if ($null -ne $data) {
        $rowCount = $data.GetLength(0)
        $columnCount = $data.GetLength(1)
    }

    $writer = [System.IO.StreamWriter]::new(
        $csvPath,
        $false,
        [System.Text.UTF8Encoding]::new($false)
    )
    try {
        $headers = for ($columnIndex = 0; $columnIndex -lt $columns.Count; $columnIndex++) {
            Convert-ToCsvField ("col_{0}_{1}" -f ($columnIndex + 1), $columns[$columnIndex].short_name)
        }
        $writer.WriteLine(($headers -join ','))

        if ($null -ne $data) {
            for ($row = $data.GetLowerBound(0); $row -le $data.GetUpperBound(0); $row++) {
                $fields = for ($column = $data.GetLowerBound(1); $column -le $data.GetUpperBound(1); $column++) {
                    Convert-ToCsvField $data[$row, $column]
                }
                $writer.WriteLine(($fields -join ','))
            }
        }
    }
    finally {
        $writer.Dispose()
    }

    $metadata = [ordered]@{
        sheet_range         = $SheetRange
        workbook            = $bookName
        worksheet           = $sheetName
        worksheet_long_name = [string]$worksheet.LongName
        origin_rows         = [int]$worksheet.Rows
        origin_columns      = [int]$worksheet.Cols
        exported_rows       = $rowCount
        exported_columns    = $columnCount
        csv_file            = [System.IO.Path]::GetFileName($csvPath)
        columns             = $columns
    }
    $metadata | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $metadataPath -Encoding utf8NoBOM
    return $metadata
}

$resolvedOutput = [System.IO.Path]::GetFullPath($OutputRoot)
$null = New-Item -ItemType Directory -Path $resolvedOutput -Force
$manifest = @()

$origin = $null
foreach ($source in $InputPath) {
    $resolvedSource = (Resolve-Path -LiteralPath $source).Path
    $projectStem = Convert-ToSafeName ([System.IO.Path]::GetFileNameWithoutExtension($resolvedSource))
    $projectOutput = Join-Path $resolvedOutput $projectStem
    $null = New-Item -ItemType Directory -Path $projectOutput -Force

    $sourceHash = (Get-FileHash -LiteralPath $resolvedSource -Algorithm SHA256).Hash.ToLowerInvariant()
    $sheetMetadata = @()
    $origin = New-Object -ComObject Origin.Application
    try {
        $origin.Visible = 0
        if (-not $origin.Load($resolvedSource)) {
            throw "Origin failed to load: $resolvedSource"
        }
        $null = $origin.Run()

        $sheetNames = Get-OriginSheetNames $origin
        foreach ($sheetName in $sheetNames) {
            $sheetMetadata += Export-OriginWorksheet $origin $sheetName $projectOutput
        }
    }
    finally {
        try { $null = $origin.Execute('doc -s') } catch { }
        try { $null = $origin.Exit() } catch { }
        [System.Runtime.InteropServices.Marshal]::FinalReleaseComObject($origin) | Out-Null
        $origin = $null
        [gc]::Collect()
        [gc]::WaitForPendingFinalizers()
    }

    $sourceHashAfter = (Get-FileHash -LiteralPath $resolvedSource -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($sourceHashAfter -ne $sourceHash) {
        throw "Source file changed during extraction: $resolvedSource"
    }

    $projectManifest = [ordered]@{
        source_file          = $resolvedSource
        source_sha256_before = $sourceHash
        source_sha256_after  = $sourceHashAfter
        extracted_at         = (Get-Date).ToString('o')
        read_only_guarantee  = 'No Save call was issued; pre/post SHA256 hashes match.'
        worksheet_count      = $sheetMetadata.Count
        worksheets           = $sheetMetadata
    }
    $projectManifestPath = Join-Path $projectOutput 'manifest.json'
    $projectManifest | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $projectManifestPath -Encoding utf8NoBOM
    $manifest += $projectManifest
}

$manifestPath = Join-Path $resolvedOutput 'manifest.json'
$manifest | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $manifestPath -Encoding utf8NoBOM
Write-Output $manifestPath
