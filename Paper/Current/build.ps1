param()
$ErrorActionPreference = 'Stop'
$paperRoot = $PSScriptRoot
$paperBuild = Join-Path $paperRoot 'build'
$paperJob = 'current-paper'
New-Item -ItemType Directory -Path $paperBuild -Force | Out-Null
Push-Location -LiteralPath $paperRoot
try {
    & pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build "-jobname=$paperJob" main.tex
    if ($LASTEXITCODE -ne 0) { throw "First LaTeX pass failed; see build/$paperJob.log." }
    & bibtex "build/$paperJob"
    if ($LASTEXITCODE -ne 0) { throw "Bibliography build failed; see build/$paperJob.blg." }
    1..2 | ForEach-Object {
        & pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build "-jobname=$paperJob" main.tex
        if ($LASTEXITCODE -ne 0) { throw "LaTeX pass failed; see build/$paperJob.log." }
    }
    $paperLog = Get-Content -LiteralPath (Join-Path $paperBuild "$paperJob.log") -Raw
    if ($paperLog -match 'There were undefined references|Citation .+ undefined|Reference .+ undefined|multiply-defined labels|Overfull') {
        throw "The PDF compiled but needs reference or layout attention; see build/$paperJob.log."
    }
    Copy-Item -LiteralPath (Join-Path $paperBuild "$paperJob.pdf") -Destination (Join-Path $paperRoot 'main.pdf')
    Copy-Item -LiteralPath (Join-Path $paperBuild "$paperJob.bbl") -Destination (Join-Path $paperRoot 'main.bbl')
    Write-Output 'Updated Current/main.pdf from Current/main.tex.'
} finally {
    Pop-Location
}
