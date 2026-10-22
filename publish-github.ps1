$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot

function Invoke-Git {
    & git @args
    if ($LASTEXITCODE -ne 0) { throw "Git failed: $args" }
}

if ((& git branch --show-current) -ne 'master') { throw 'Run this script on the master branch.' }
if (& git ls-files --unmerged) { throw 'Resolve the existing merge conflicts before publishing.' }

if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
    throw 'Install GitHub CLI first: winget install --id GitHub.cli --exact . Then reopen PowerShell and run this script again.'
}
& gh auth status
if ($LASTEXITCODE -ne 0) {
    & gh auth login --hostname github.com --git-protocol https --web
    if ($LASTEXITCODE -ne 0) { throw 'GitHub login failed.' }
}
& gh auth setup-git
if ($LASTEXITCODE -ne 0) { throw 'Git authentication setup failed.' }

$accountJson = & gh api user
if ($LASTEXITCODE -ne 0) { throw 'Cannot read GitHub account.' }
$account = $accountJson | ConvertFrom-Json
if ($account.login -ne 'Mavis-Mai') { throw 'Please log in to the Mavis-Mai GitHub account.' }
if (-not (& git config user.name)) { Invoke-Git config user.name $account.login }
if (-not (& git config user.email)) {
    Invoke-Git config user.email "$($account.id)+$($account.login)@users.noreply.github.com"
}

Invoke-Git lfs install --local
Invoke-Git rm -r --cached --ignore-unmatch -- ':(glob)**/__pycache__/**' ':(glob)**/*.pyc'
Invoke-Git add --all
& git diff --cached --quiet
if ($LASTEXITCODE -eq 1) {
    Invoke-Git commit -m 'Rename project to GSAeval and update code and datasets'
} elseif ($LASTEXITCODE -ne 0) { throw 'Cannot inspect staged changes.' }

# Resolve by repository ID after renaming so a pre-existing GSAeval repository is never reused.
$oldJson = & gh api repos/Mavis-Mai/GSAqualifier 2>$null
if ($LASTEXITCODE -eq 0) {
    $oldRepo = $oldJson | ConvertFrom-Json
    $targetBranch = $oldRepo.default_branch
    if ($oldRepo.name -eq 'GSAqualifier') {
        & gh api --method PATCH "repositories/$($oldRepo.id)" -f name=GSAeval
        if ($LASTEXITCODE -ne 0) { throw 'Repository rename failed. No code was pushed.' }
    } elseif ($oldRepo.name -ne 'GSAeval') { throw 'Unexpected repository name.' }
} else {
    throw 'Cannot resolve the original repository. Check access to Mavis-Mai/GSAqualifier.'
}
Invoke-Git remote set-url origin https://github.com/Mavis-Mai/GSAeval.git
Invoke-Git fetch origin
& git show-ref --verify --quiet "refs/remotes/origin/$targetBranch"
if ($LASTEXITCODE -eq 0) {
    Invoke-Git merge --no-edit --allow-unrelated-histories "origin/$targetBranch"
}
Invoke-Git push --set-upstream origin "HEAD:refs/heads/$targetBranch"
if ((& git rev-parse HEAD) -ne (& git rev-parse "origin/$targetBranch")) {
    throw 'Remote verification failed.'
}
Write-Host 'Published successfully: https://github.com/Mavis-Mai/GSAeval'
