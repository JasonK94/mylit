# Project-Specific Context

## Primary Goal
Deliver a reusable R toolkit (`myR`) for single-cell RNA-seq and GeoMx spatial transcriptomics workflows.  
The package should bundle proven analysis modules—pseudobulk DEG, pathway enrichment, pseudotime, cell–cell interaction, GeoMx preprocessing—and expose them through documented, reproducible interfaces.

## AI Model
Anything (focus on code reasoning, documentation drafting, and data-processing support).

---

# General Context for AI Assistant

## Project Initialization Workflow

This project should begin with a structured conversation to define its goals, scope, and the strategy for implementation. Follow these steps:

1.  **User's Goal Statement**: Gather requirements tied to `myR` features or documentation gaps.
2.  **AI-led Scoping and Strategy Discussion**:
    *   Clarify which module or document (e.g., pseudobulk, Milo, README) needs attention.
    *   Identify relevant commits, scripts, or vignettes (`docs/functions/function_analysis.md`, `R/*.R`).
    *   Agree on deliverables (code changes, docs, tests).
3.  **Create a Detailed Plan**: Produce a TODO list before editing, referencing target files.
4.  **Begin Implementation**: Execute plan incrementally, verifying results at each step.

## Evolving This Context File

This document is not static. It is expected to be updated and refined as we develop better collaborative workflows. The process for updating it is as follows:

1.  **Propose a Change**: Either the user or the AI can propose a change to this document based on a new idea for improving the process.
2.  **Document Failures**: If the AI makes a critical mistake (like overwriting a file or using a deprecated value), the user or AI must add a new rule to the "Project-Specific Caveats" section to prevent it from happening again.
3.  **Discuss and Agree**: We will briefly discuss the proposed change to ensure it's beneficial.
4.  **Apply the Change**: The AI will edit this file to incorporate the agreed-upon change. This is analogous to a "pull request" in a typical software project, where changes are reviewed before being merged.

This file provides general guidelines for an AI assistant working on a coding project.

## Core Principles
1.  **Understand the Goal**: Tie every change back to the package roadmap (e.g., Milo enhancements, documentation sync).
2.  **Plan Your Work**: Create and maintain TODO lists; group related doc/code updates.
3.  **Be Systematic**: Update one module or document at a time; keep commits focused.
4.  **Explain Your Actions**: Describe why commands or edits are required (e.g., “inspect `b260bb1` for Milo context”).
5.  **Self-Correction**: If a command fails or findings conflict, summarize the issue and propose a revised plan.

## 🪲 Project-Specific Caveats (Learned Lessons)

*This section is a living document that records project-specific rules and lessons learned from past mistakes. Both the AI and the user are responsible for updating this list to prevent repeated errors.*

---
*(Add project-specific rules here as they are discovered.)*
- Milo 파이프라인은 `nhoods(milo)` 희소 행렬을 요약할 때 계산 비용이 크다. 가능하면 `save=TRUE`로 캐시를 활용하고, 강제 재계산이 필요한 경우에만 `force_run`을 사용한다.
- `plotNhoodGraphDA()` 실행 전에는 반드시 `buildNhoodGraph(milo)`가 수행돼야 한다. 캐시에서 읽어온 Milo 객체의 그래프 슬롯은 비어 있을 수 있으므로 재생성한다.
- SpatialFDR이 모두 `alpha` 이상일 때 `plotDAbeeswarm()`이 실패한다. 기본 fallback 메트릭(`PValue` 등)을 준비하고, 색상 구간을 rank 기반으로 재설정할 수 있도록 한다.
- 모든 Milo 관련 R 세션은 `st/start.R`에서 `renv`를 활성화한 뒤 시작한다. 루트 경로에서 바로 R을 실행하면 필수 패키지가 로드되지 않아 디버깅이 불가능해진다.
- `run_milo_pipeline()` 같은 핵심 함수는 반드시 `myR/R/` 최상위에 두고, 하위 서브폴더(예: `myR/R/analysis/`)에 넣지 않는다. 그렇지 않으면 `devtools::load_all()`이 함수를 namespace에 노출하지 못한다.
- `SingleCellExperiment::colData<-`처럼 네임스페이스를 강제하는 setter는 사용하지 않는다. 특정 버전에서 export되지 않아 `'not an exported object'` 오류가 발생하므로 S4Vectors/SummarizedExperiment accessor만 사용한다.
- GEM barcode suffix(`-1`, `-2`)를 제거하거나 샘플 ID를 자의적으로 수정하지 않는다. 과거 AI가 suffix를 제거했다가 sample ID가 중복되어 DA 전 단계가 실패했다.
- 자동화된 `git switch`는 사용하지 않는다. 다른 worktree를 잘못 가리켜 작업물이 사라질 뻔한 사례가 있다.
- 캐시 `.qs`를 재사용할 때는 동일 데이터/파라미터인지 먼저 확인한다. 필요하면 `cache_files` 인자를 통해 명시적으로 경로를 지정하고, `milo$commands` 로그를 확인해 잘못된 캐시를 재활용하지 않도록 한다.

## Code Quality & Style
1.  **Clean & Readable**: Write clean, well-structured, and commented code. Follow the existing coding style of the project.
2.  **Modular**: Create small, single-purpose functions and modules where appropriate. Avoid monolithic scripts.
3.  **Configuration over Hardcoding**: Use configuration files (`.json`, `.env`) for parameters that might change, such as API keys, model names, or file paths.
4.  **Error Handling**: Implement robust error handling. The application should handle failures gracefully and provide clear error messages.

## Project Management
1.  **Version Control**: Keep documentation changes (DEVLOG/CHANGELOG) synchronized with relevant commits.
2.  **Documentation**: `README.md` (English) and `README_Korean.md` describe capabilities; update them when modules evolve.
3.  **Dependency Management**: Reflect dependency changes in `DESCRIPTION` and note rationale in CHANGELOG.

## Communication
1.  **Clarity and Conciseness**: Be clear and to the point. Avoid jargon where possible.
2.  **Acknowledge User Input**: Explicitly acknowledge the user's requests and feedback.
3.  **Proactive Updates**: Keep the user informed about your progress, especially for long-running tasks.

## Development History & Push Policy

- Maintain two documents at repo root:
  - `DEVLOG.md`: narrative context of sessions, decisions, and next steps.
  - `CHANGELOG.md`: semantic, versioned record of notable changes.
- Workflow:
  1. Record intent and context in `DEVLOG.md` before significant work.
  2. Make local commits frequently; push only major/stable changes or when collaboration requires it.
  3. Summarize major changes in `CHANGELOG.md` using Keep a Changelog style.
  4. Reference the relevant DEVLOG entry in commit messages when helpful.

## Project Artifacts

A new project initialized via `cinit` contains several key files. Understand their roles:

1.  **`context.md`**: The **Single Source of Truth** for the AI agent. This is your primary tool for guiding the AI. It defines the project's goal, scope, and technical constraints. It should be updated continuously as the project evolves.
2.  **`NEXT_STEPS.md`**: A **bootstrapping guide**. It contains the ideal first prompt to give the AI agent to kickstart the development process in a structured way. Its purpose is fulfilled after this first prompt.
3.  **`DEVLOG.md`**: The **project's narrative log**. Use it to record the "why" behind decisions, track experiments, and maintain context between development sessions. This is for both you and the AI to review.
4.  **`CHANGELOG.md`**: The **formal record of changes**. Use it to document notable updates, bug fixes, and new features, typically adhering to Semantic Versioning.
