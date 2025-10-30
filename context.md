# Project-Specific Context

## Primary Goal
the data is acquired from UC patient's endoscopic biopsy sample of active resion in rectum using GeoMx DSP platform. patients are on biologics treatment (Infliximab, Ustekinumab, Vedolizumab). samples are obtained pre/post treatment. Most of but not all samples have one ROI which has two AOIs(PanCK positive / negative). Most patients got just one drug treatment, but some tried two on different periods. The goainvestigate each drug's molecular effect, the difference between each drug, and the difference between clinical responders and nonresponders. All the aanalysis should be done for each tissue type(AOI type; PanCK positive / negative; epithelium / stroma). patient ID: $emrid, sample ID: $patho_id, ROI id:

## AI Model
anything.

---

# General Context for AI Assistant

## Project Initialization Workflow

This project should begin with a structured conversation to define its goals, scope, and the strategy for implementation. Follow these steps:

1.  **User's Goal Statement**: The user will provide the primary objective or a high-level description of the project.
2.  **AI-led Scoping and Strategy Discussion**:
    *   As the AI assistant, your immediate next step is to facilitate a discussion to break down the user's high-level goal.
    *   Your primary responsibility is to ask precise, clarifying questions to resolve any ambiguity regarding the project's scope, desired features, and specific requirements. Do not proceed if the goal is not well-defined.
    *   Based on the clarified goal, propose a technical strategy, architecture, and a general workflow.
    *   Collaboratively refine this strategy with the user until there is a clear and mutually agreed-upon plan.
3.  **Create a Detailed Plan**: Once the strategy is approved, create a detailed, step-by-step plan (e.g., a TODO list) that outlines the tasks required for execution.
4.  **Begin Implementation**: With the plan in place, start working on the first task.

## Evolving This Context File

This document is not static. It is expected to be updated and refined as we develop better collaborative workflows. The process for updating it is as follows:

1.  **Propose a Change**: Either the user or the AI can propose a change to this document based on a new idea for improving the process.
2.  **Document Failures**: If the AI makes a critical mistake (like overwriting a file or using a deprecated value), the user or AI must add a new rule to the "Project-Specific Caveats" section to prevent it from happening again.
3.  **Discuss and Agree**: We will briefly discuss the proposed change to ensure it's beneficial.
4.  **Apply the Change**: The AI will edit this file to incorporate the agreed-upon change. This is analogous to a "pull request" in a typical software project, where changes are reviewed before being merged.

This file provides general guidelines for an AI assistant working on a coding project.

## Core Principles
1.  **Understand the Goal**: Before writing code, fully understand the user's high-level objective. Ask clarifying questions if the goal is ambiguous.
2.  **Plan Your Work**: For any non-trivial request, create a plan and share it. Use a TODO list to track progress and mark items as complete.
3.  **Be Systematic**: Make one logical change at a time. Verify each change before moving to the next. Avoid making many unrelated changes in a single step.
4.  **Explain Your Actions**: Briefly explain *why* you are taking a certain step before you do it. Provide the commands for the user to run, explaining what each one does.
5.  **Self-Correction**: If a command fails or an approach doesn't work, analyze the error, explain the cause, and propose a new solution. Don't repeat the same mistake.

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
1.  **Version Control**: Use Git for version control. Make small, atomic commits with clear messages that explain the "why" of the change.
2.  **Documentation**: Keep `README.md` and other documentation up-to-date. Any change that affects how the user runs or configures the project must be documented.
3.  **Dependency Management**: Use a package manager (`package.json`, `requirements.txt`, etc.) and keep dependencies clean. Explain why a new dependency is needed before adding it.

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
