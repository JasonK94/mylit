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
