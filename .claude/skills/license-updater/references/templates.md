# License File Templates

Use these templates when the generation script is unavailable or the situation
requires custom formatting.

## THIRD_PARTY_LICENSES

```
Third-Party Licenses
====================

This file documents the licenses for third-party code included in or adapted
by this project.


================================================================================
{PROJECT_NAME}
================================================================================

Source:       {PROJECT_URL}
License:      {SPDX_IDENTIFIER}
Files:        {PORTED_FILES_PATH}

{DESCRIPTION_OF_WHAT_WAS_PORTED}

{PROVENANCE_CHAIN — if the code has a history of forks/rewrites, list each
step with author and license}

{FULL_ORIGINAL_LICENSE_TEXT — copy verbatim, do not paraphrase}
```

When appending to an existing THIRD_PARTY_LICENSES file, add a new section
starting from the `====` separator line. Do not modify existing sections.

## Header template (for TypeScript / JavaScript files)

```typescript
/**
 * {BRIEF_DESCRIPTION — one line saying what this code does}
 *
 * {PORT_DESCRIPTION — e.g. "TypeScript port of {project} ({url})"}
 * {ORIGINAL_ALGORITHM_CREDIT — if applicable}
 *
 * {ORIGINAL_COPYRIGHT_LINES — one per line, all of them}
 * {TARGET_COPYRIGHT_LINE}
 *
 * Licensed under the {SPDX_IDENTIFIER} License. See THIRD_PARTY_LICENSES
 * for the full license text and provenance chain of the original code.
 */
```

### Header guidelines

- Keep under 15 lines total
- Include ALL original copyright holders
- Reference THIRD_PARTY_LICENSES instead of embedding the full license
- For BSD-3-Clause, do NOT include the full "neither the name..." clause
  in the header — the THIRD_PARTY_LICENSES file covers this

## Header template (for Python files)

```python
# {BRIEF_DESCRIPTION}
#
# {PORT_DESCRIPTION}
# {ORIGINAL_ALGORITHM_CREDIT}
#
# {ORIGINAL_COPYRIGHT_LINES}
# {TARGET_COPYRIGHT_LINE}
#
# Licensed under the {SPDX_IDENTIFIER} License. See THIRD_PARTY_LICENSES
# for the full license text and provenance chain.
```

## Header template (for C / C++ files)

```c
/*
 * {BRIEF_DESCRIPTION}
 *
 * {PORT_DESCRIPTION}
 * {ORIGINAL_ALGORITHM_CREDIT}
 *
 * {ORIGINAL_COPYRIGHT_LINES}
 * {TARGET_COPYRIGHT_LINE}
 *
 * Licensed under the {SPDX_IDENTIFIER} License. See THIRD_PARTY_LICENSES
 * for the full license text and provenance chain.
 */
```

## README License section

Add the following beneath the existing license statement in the project's
`README.md`. If the README has no `## License` section yet, create one.

```markdown
## License

This project is licensed under the [{TARGET_LICENSE}]({LICENSE_FILE_PATH}).

{PORTED_COMPONENT_DESCRIPTION} is a {LANGUAGE} port of
[{PROJECT_NAME}]({PROJECT_URL}).
See [THIRD_PARTY_LICENSES](./THIRD_PARTY_LICENSES) for details on third-party
code provenance.
```

### README License guidelines

- Keep it to 2–3 sentences maximum
- Always link to both the main `LICENSE` file and `THIRD_PARTY_LICENSES`
- Name the specific component that was ported (e.g., "The LSODA solver")
- Name the upstream project and link to it
- If multiple third-party components exist, list each on its own line

## LICENSE_GUIDE.md

```markdown
# License Attribution Guide

## Overview

This guide explains how to properly attribute the third-party code ported
into this project from [{PROJECT_NAME}]({PROJECT_URL}).

## Files to add/modify

### 1. `THIRD_PARTY_LICENSES` → repo root

{Explain: new file or append to existing}

### 2. Header in each ported file → `{PORTED_FILES_PATH}/*.{ext}`

{Show the header template}

### 3. `LICENSE` → no changes needed

{Explain why the main license stays as-is}

### 4. `README` → optional update

{Suggest a license section addition if appropriate}

## License compatibility

Upstream: {UPSTREAM_LICENSE}
Target:   {TARGET_LICENSE}
Status:   ✅ Compatible / ❌ Incompatible

{Brief explanation of compatibility and any special conditions}
```
