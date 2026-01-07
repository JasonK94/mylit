#!/bin/bash
# Documentation subtree synchronization script
# Usage: ./scripts/docs-sync.sh [push|pull] [branch]

set -e

ACTION=${1:-pull}
BRANCH=${2:-main}
REPO_ROOT="/home/user3/data_user3/git_repo/mylit"

# Change to repo root (works from any subdirectory)
cd "$REPO_ROOT"

# Check if docs-main-origin remote exists
if ! git remote | grep -q "^docs-main-origin$"; then
    echo "Error: docs-main-origin remote not found"
    echo "Please run setup script first: ./scripts/setup-docs-subtree.sh <repo-url>"
    exit 1
fi

if [ "$ACTION" = "push" ]; then
    echo "Pushing docs-main changes to remote ($BRANCH)..."
    git subtree push --prefix=docs-main docs-main-origin $BRANCH
    echo "✓ Successfully pushed docs-main to remote"
    
elif [ "$ACTION" = "pull" ]; then
    echo "Pulling latest docs-main from remote ($BRANCH)..."
    git subtree pull --prefix=docs-main docs-main-origin $BRANCH
    echo "✓ Successfully pulled latest docs-main from remote"
    
else
    echo "Usage: $0 [push|pull] [branch]"
    echo ""
    echo "Examples:"
    echo "  $0 pull main    # Pull latest docs from remote"
    echo "  $0 push main    # Push local changes to remote"
    exit 1
fi

