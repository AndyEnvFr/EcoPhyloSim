# Contributing to EcoPhyloSim

Contributions to EcoPhyloSim and its R package PhyloSim are welcome! This document describes how to contribute effectively and what to expect from the process. Tasks of the EcoPhyloSim team (maintainers and collaborators with write access to this repository) are described in [Maintaining EcoPhyloSim](MAINTAINING.md).

## Ways to contribute

### Asking for help and reporting issues

PhyloSim is open-source software, and in the spirit of open source, questions and answers should be public, so that other people can profit from them as well. Please ask questions, report problems and suggest features on the [EcoPhyloSim issue tracker](https://github.com/TheoreticalEcology/EcoPhyloSim/issues), after searching the existing issues. If you have concerns posting to a public forum, you can also write to the maintainer by e-mail (run `maintainer("PhyloSim")` in R to get the address).

It is easier to help you if your question is well structured. Typically, it contains:

1. the motivation of the study,
2. your PhyloSim code as you execute it (e.g., the parameters passed to `createCompletePar()`),
3. PhyloSim results / error messages, and the output of `sessionInfo()`,
4. your question.

For problems and errors, please always include a minimal reproducible example: code that runs as is in a fresh R session and shows the problem. See [How to make a great R reproducible example](https://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example) and the [reprex](https://cran.r-project.org/package=reprex) package. Keep simulations small and fast (small grid, few generations) and set a seed. Problem reports without a reproducible example may be closed.

Please write issues so that they are understandable without earlier discussions. If you have already diagnosed the problem or have ideas for a fix, including that analysis is very helpful. You are also welcome to comment on existing issues, for example by adding examples or proposing solutions.

### Submitting pull requests

Please open or comment on an issue before you start working on a pull request, so that the approach can be discussed before you invest time in an implementation. This avoids situations where a pull request cannot be merged because of design conflicts that could have been resolved earlier. This applies in particular to changes of the simulation model (C++ code in `PhyloSim/src/`), which change the results of existing simulations. Only trivial fixes (e.g., typos or a factual error in the documentation) can go directly to a pull request. See [Pull request workflow](#pull-request-workflow) below.

## Policy on AI-generated contributions

This policy applies to contributors outside the EcoPhyloSim team; only the disclosure of AI use (last paragraph) applies to everyone.

**Issues, pull requests and comments created autonomously by AI agents or automated tools will not be accepted.** They will be closed or deleted without review.

The maintainer must be able to understand, evaluate, and take long-term responsibility for every change integrated into PhyloSim. AI-generated code frequently appears correct but introduces subtle errors that are difficult to detect without deep familiarity with the codebase, and the burden of reviewing such contributions falls disproportionately on the maintainer. For a scientific simulation model where correctness is critical, this risk is not acceptable.

This policy covers issues, pull requests and comments submitted autonomously by AI agents, as well as contributions where the submitter does not have a full, independent understanding of the proposed changes. Having a tool post a contribution at your explicit request, after you have checked its content, does not count as autonomous.

Using AI tools as a personal aid during development is a matter of individual choice. However, any contribution must be authored and understood by a human who can explain, defend, and iterate on the proposed changes.

If you are unsure whether your contribution complies with this policy, ask yourself: can you explain every line of your proposed change and its consequences for the rest of the codebase? If not, please do not submit it as a pull request — an issue describing the problem or idea is a better starting point.

For everyone, including the EcoPhyloSim team: if AI tools wrote substantial parts of your code or text, please mention this briefly in your issue, pull request, comment or commit (one sentence is enough, e.g., "Code drafted with an AI assistant, checked and tested by me."; in commits, use a trailer such as `Assisted-by: <tool> (<model>)`). Instructions for AI coding agents, including those used by the EcoPhyloSim team, are in [AGENTS.md](AGENTS.md).

## Pull request workflow

`master` is the development branch; releases are tagged on `master` (e.g., `v0.3`). There are no separate development or release branches.

- Members of the EcoPhyloSim team create branches in this repository; other contributors work in a fork, where the branch name is up to them.
- Every branch in this repository addresses at least one issue (except trivial fixes, see [Submitting pull requests](#submitting-pull-requests)). Create it from `master` and name it as follows, with `-` as separator and `<shortName>` in lowerCamelCase:
  - one issue: `i<issue number>-<shortName>`, e.g., `i115-sacLogSpacing`.
  - several closely related issues for the same release: `v<release>-<shortName>`, e.g., `v0.4-docFixes`.
  - trivial fixes without an issue: `<shortName>`, e.g., `fixTypoVignette`.
  - Do not name a branch exactly like a release tag (e.g., `v0.4`).
- Keep each pull request focused on one issue, or on several closely related issues for the same release, and open it against `master`. Do not split one issue across several pull requests; follow-up work gets a new issue.
- Keep the diff small: no unrelated reformatting (indentation, line endings, whitespace) or refactoring, and reuse existing helper functions instead of duplicating logic. If reformatting is needed, it goes into a separate pull request without other changes.
- Never commit compiled files (`.o`, `.so`, `.dll`) or other build output, and do not remove the entries for them in `.gitignore`.
- Once the approach is agreed in the issue, open the pull request early, as a draft if the work is not finished yet, so that work in progress is visible and linked to its issues. If the work is abandoned, close the pull request with a short explanation.
- A pull request or commit that resolves an issue should contain `Fixes #115` or `Closes #115` (with your issue number) in its description or commit message, so that GitHub closes the issue automatically when it is merged into `master`. For several issues, repeat the keyword for each, e.g., `Fixes #115, fixes #116`.
- The pull request description has the following sections; keep it short, details belong in the issues:

  ```markdown
  **Problem addressed:** one or two sentences, with references to the issues

  **Changes**
  - one short bullet per change, with "Fixes #n" for each issue it resolves

  **Open problems**
  - list (or "none")

  **Checks and tests**
  - what was checked and which tests were added (and whether they fail without the change); what was not checked
  ```

  Pull requests from contributors outside the EcoPhyloSim team are reviewed before they are merged.
- Titles of issues and commit subject lines are short plain text (ideally under 70 characters, no backticks) naming the function or area and the problem or change, so that a list of titles is understandable without opening them, e.g., "sac(): default plot sizes are not log-spaced"; avoid vague titles such as "Small fixes" or "Revisions to X".
- In commit messages, separate the body (if any) from the subject line by an empty line.
- Pull request titles are short plain text stating the general idea, without issue numbers, e.g., "Fixes to grouping calculations".
- After the merge, branches in this repository are deleted automatically.

## Code and documentation standards

PhyloSim is an R package with a simulation model in C++ (via Rcpp); the package itself is in the subfolder `PhyloSim/`. Please follow these conventions for new and changed code (some existing code does not yet follow them):

- Use meaningful names; function and variable names in lowerCamelCase. Accessor functions start with `get` (e.g., `getSpecies()`, `getSpeciesTraits()`).
- Help files and `NEWS.md` are written in American English; parameter descriptions start in lower case.
- Code that addresses a specific problem described in an issue and is not self-explanatory gets a short comment with the issue number and, if it fits in a few words, the problem it solves, e.g., `# fixes #115 (log-spaced plot sizes)`. Use such comments sparingly, typically in one or two places per issue.

### Documentation and examples

- New or changed functions should include roxygen2 documentation. The help files in `PhyloSim/man/` and `PhyloSim/NAMESPACE` are generated from it; do not edit them by hand, but run `devtools::document("PhyloSim")` from the repository root, and commit only the changes that belong to your change.
- Examples always go into a separate file in `PhyloSim/inst/examples/`, named after the R file that contains the function, and are included in the roxygen block with `@example`, e.g., `#' @example inst/examples/plotSAC-help.R` for `sac()` in `PhyloSim/R/plotSAC.R`. Do not write examples directly in the roxygen block with `@examples`. Examples must run in a few seconds; wrap longer ones in `\donttest{}`.

### C++ code

- The simulation model is in `PhyloSim/src/`. Print with `Rcpp::Rcout` (and `Rcpp::Rcerr`), never with `std::cout`, `std::cerr` or `printf`; signal errors with `Rcpp::stop()`, never with `exit()` or `abort()`.
- `PhyloSim/src/RcppExports.cpp` and `PhyloSim/R/RcppExports.R` are generated: after changing an `// [[Rcpp::export]]` function, run `Rcpp::compileAttributes("PhyloSim")`. `PhyloSim/src/RcppExportsIni.c` registers `callModel` and is maintained by hand: if the number of arguments of `callModel` changes, update both the `SEXP` declaration and the count in `CallEntries`.
- `pugixml.cpp`, `pugixml.hpp` and `pugiconfig.hpp` are bundled third-party code; do not modify them.
- Changes to the model that change simulation results, even for the default parameters, must be stated in the pull request description (and, as every user-visible change, in `NEWS.md`).

### Tests and checks

- If your contribution changes or adds functionality, please include or update tests in `PhyloSim/tests/testthat/`. Tests should check behavior that matters, not cosmetic details, and should be fast; mark slow tests with `skip_on_cran()`. Model changes should be tested with a fixed seed.
- Tests, examples and the vignette use at most 2 cores (e.g., `runSimulationBatch(..., parallel = 2)`), never `parallel = TRUE` or `"auto"`, as required by CRAN.
- A CRAN release is planned, so `R CMD check --as-cran` should be clean; a pull request must not add ERRORs, WARNINGs or NOTEs. `R CMD check` runs on pull requests to `master` (GitHub Actions; for first-time contributors after approval by the EcoPhyloSim team) and must pass. Locally, it is usually sufficient to run the tests related to your change, e.g., `testthat::test_file("PhyloSim/tests/testthat/test-runSimulation.R")` after `devtools::load_all("PhyloSim")`.

### Dependencies

- Discuss new dependencies in an issue first. Prefer base R or packages already in `PhyloSim/DESCRIPTION`; packages used only in optional functions go into `Suggests` and are checked with `requireNamespace()`. Do not call `library()` or `require()` in package code.

### NEWS

- `PhyloSim/NEWS.md` must always be updated: every pull request or commit with a user-visible change to the package adds, in the same pull request or commit, a short (ideally one-line) entry under the matching subheading (New features, Changes, Bug fixes, Documentation) of the topmost version heading (no entry for bugs introduced after the last release), ending with the issue or pull request number(s) if there are any, e.g., ``* Default plot sizes in `sac()` are log-spaced, as documented (#115).`` For substantial contributions from outside the EcoPhyloSim team (see [Attribution and authorship](#attribution-and-authorship)), add the author's GitHub username first, e.g., `(@user, #115)`.

## Attribution and authorship

Contributors who make substantial contributions to PhyloSim will be acknowledged by name in `PhyloSim/NEWS.md` and/or as contributors in `PhyloSim/DESCRIPTION`. Substantial contributions include, but are not limited to: adding a new model mechanism or a meaningful new function or test, identifying and fixing a non-trivial bug, or making significant improvements to the documentation or vignette.

Minor corrections such as fixing typos, correcting punctuation, or making small formatting adjustments do not typically result in named acknowledgement, as is standard in open-source scientific software. If you are uncertain whether a planned contribution would qualify, you are welcome to ask in the relevant issue.
