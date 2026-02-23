# Contributing to diff-grok 🚀

Thank you for your interest in contributing to diff-grok!  
We welcome bug reports, documentation improvements, new features, and help with tests or examples.

## ✅ How to contribute

### Reporting issues

- Please search existing issues before creating a new one — your problem or idea may already be discussed or solved.  
- When reporting a bug, include as much information as possible: error messages, minimal reproducible example, environment (Node version, OS, TypeScript version), and expected vs actual behavior.  
- For feature requests — describe the motivation, proposed API (if relevant), and possible use-cases.

### Submitting pull requests (PR)

1. Fork the repository and create a branch:  
   `feat/your-feature` or `fix/bug-description` or `docs/update-readme`  
2. Make sure your code follows the existing style (TypeScript, ESLint, Prettier — as used in the repo).  
3. Add or update tests/examples if applicable.  
4. Run linter and tests before committing:

```bash
   npm install  
   npm run lint  
   npm test
```

5. Commit with a clear message. Prefer descriptive messages, e.g.:

```bash
feat: add new Rosenbrock method,
fix: correct edge behaviour in solver,
docs: improve README example for ODE usage.
```

6. Rebase or merge main upstream if needed; make sure no conflicts remain.

7. Create a PR, describing what you’ve changed and referencing related issues (if any).

## 🧑‍💻 Code style & conventions

Project uses TypeScript.

Follow the existing coding conventions and folder structure.

Keep changes focused and minimal — avoid mixing unrelated changes (e.g. code + formatting + docs) in one PR.

Write clear, self-contained commits.

## 📚 Documentation & Examples

If you fix or add public API — update README or add an example.

For new features — prefer adding minimal examples or tests demonstrating usage.

If you change behavior in edge cases — mention it clearly in docs/README.

## 🗣 Communication & Etiquette

Use GitHub Issues and Pull Requests for all discussions.

Be respectful and patient: contributors and maintainers may be in different timezones or busy.

Before opening an Issue or PR — read existing content (README, issues, pull requests) to avoid duplicates.

## 🔧 CI, Tests, and Quality

Ensure all existing tests pass; add new tests if adding functionality.

Code should pass linting/formatting checks.

If your change is substantial — consider adding a test or example illustrating it.

Thank you for helping make diff-grok better! 🎉
