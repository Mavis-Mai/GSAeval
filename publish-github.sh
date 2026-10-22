#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if ! command -v gh >/dev/null && [[ -x /tmp/gsaeval-gh/bin/gh ]]; then
    export PATH="/tmp/gsaeval-gh/bin:$PATH"
fi
command -v gh >/dev/null || { echo 'Install GitHub CLI (gh) in WSL first.' >&2; exit 1; }
[[ $(git branch --show-current) == master ]] || { echo 'Expected master branch.' >&2; exit 1; }
[[ -z $(git ls-files --unmerged) ]] || { echo 'Resolve existing merge conflicts first.' >&2; exit 1; }
gh auth status || gh auth login --hostname github.com --git-protocol https --web
gh auth setup-git
login=$(gh api user --jq .login)
[[ "$login" == Mavis-Mai ]] || { echo 'Please log in as Mavis-Mai.' >&2; exit 1; }
git config user.name >/dev/null || git config user.name "$login"
git config user.email >/dev/null || git config user.email "$(gh api user --jq .id)+${login}@users.noreply.github.com"
git lfs install --local
git rm -r --cached --ignore-unmatch -- ':(glob)**/__pycache__/**' ':(glob)**/*.pyc'
git add --all
git update-index --chmod=+x GSAeval
if ! git diff --cached --quiet; then
    git -c diff.renames=false commit -m 'Rename project to GSAeval and update code and datasets'
fi
repo_id=$(gh api repos/Mavis-Mai/GSAqualifier --jq .id)
repo_name=$(gh api "repositories/$repo_id" --jq .name)
target_branch=$(gh api "repositories/$repo_id" --jq .default_branch)
if [[ "$repo_name" == GSAqualifier ]]; then
    gh api --method PATCH "repositories/$repo_id" -f name=GSAeval >/dev/null
elif [[ "$repo_name" != GSAeval ]]; then
    echo 'Unexpected repository name.' >&2
    exit 1
fi
git remote set-url origin https://github.com/Mavis-Mai/GSAeval.git
git fetch origin
if git show-ref --verify --quiet "refs/remotes/origin/$target_branch"; then
    git merge --no-edit --allow-unrelated-histories "origin/$target_branch"
fi
git push --set-upstream origin "HEAD:refs/heads/$target_branch"
[[ $(git rev-parse HEAD) == $(git rev-parse "origin/$target_branch") ]]
echo 'Published successfully: https://github.com/Mavis-Mai/GSAeval'
