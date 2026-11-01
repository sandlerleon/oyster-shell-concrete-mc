# Publish the repository to GitHub (run in PowerShell; requires git and the GitHub CLI `gh`, already logged in as sandlerleon)
$repoDir = "C:\Users\Leon\Downloads\Oyster_work\repo_p1"
$slug    = "sandlerleon/oyster-shell-concrete-mc"
Set-Location $repoDir

function Check($what) { if ($LASTEXITCODE -ne 0) { Write-Host "FAILED: $what" -ForegroundColor Red; exit 1 } }

# 1. create the public repository if it does not exist yet (a 'not found' message from gh is expected the first time)
$exists = $false
$prev = $ErrorActionPreference; $ErrorActionPreference = "Continue"
$null = & gh repo view $slug 2>&1
if ($LASTEXITCODE -eq 0) { $exists = $true }
$ErrorActionPreference = $prev
if (-not $exists) {
    gh repo create $slug --public --description "Oyster-shell concrete Monte Carlo: coupled hydration-transport-damage-carbon model with out-of-sample validation (Constr. Build. Mater. manuscript)"
    Check "gh repo create"
}

# 2. remote, push branch and the tags (v1.0.0 is the archived first version, v1.1.0 the current one)
$ErrorActionPreference = "Continue"
git remote remove origin 2>$null
$ErrorActionPreference = $prev
git remote add origin "https://github.com/$slug.git"; Check "git remote add"
git branch -M main
git push -u origin main; Check "git push main"
git push origin --tags --force; Check "git push tags"

# 3. release with the manuscript PDFs
gh release create v1.1.0 `
    "manuscript\Oyster_P1_CBM_Manuscript.pdf" `
    "manuscript\Oyster_P1_CBM_Supplementary_Information.pdf" `
    "manuscript\Cover_Letter_CBM.pdf" `
    --repo $slug --title "v1.1.0" `
    --notes "Code, tests, results and manuscript for the Construction and Building Materials submission. Zenodo: software 10.5281/zenodo.23867048, preprint 10.5281/zenodo.23867044."
Check "gh release create"

Write-Host "Done: https://github.com/$slug" -ForegroundColor Green
