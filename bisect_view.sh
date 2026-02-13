#!/bin/bash
# Bisect script to find where View() broke

TEST_COMMITS=(
    "944ddab"  # 작동함 (확인됨)
    "4aa6eae"  # 다음 커밋
    "077bc10"  # merged plots-dev
    "09f86cd"  # merged milo  
    "4754f14"  # merged pipe
    "0e844ff"  # merged pipe2
    "bebfcec"  # docs: add MASC CLI usage
    "6824d59"  # printmarkerforllm moved
    "06ffd4c"  # Fix: Move S4Vectors (현재 main)
)

echo "=== Testing View() across commits ==="
echo ""

for commit in "${TEST_COMMITS[@]}"; do
    echo "----------------------------------------"
    echo "Testing commit: $commit"
    echo "----------------------------------------"
    
    git checkout "$commit" 2>&1 | grep -E "(HEAD|Already)"
    
    # Run test
    Rscript -e "
    suppressMessages({
        devtools::load_all('myR', quiet=TRUE)
        library(Seurat, quietly=TRUE)
    })
    counts <- matrix(rpois(200, 5), 10, 20)
    rownames(counts) <- paste0('G', 1:10)
    colnames(counts) <- paste0('C', 1:20)
    obj <- CreateSeuratObject(counts=counts, project='Test')
    
    result <- tryCatch({
        utils::View(obj)
        'SUCCESS'
    }, error = function(e) {
        paste('FAILED:', conditionMessage(e))
    })
    
    cat(sprintf('%s: %s\n', '$commit', result))
    " 2>&1 | tail -1
    
    echo ""
done

echo "=== Bisect complete ==="
