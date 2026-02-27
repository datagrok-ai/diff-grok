---
name: license-updater
description: >
  Update and manage third-party license attribution when porting code from
  external libraries. Use this skill whenever the user mentions porting,
  translating, or adapting code from another library (e.g. C/C++/Python to
  TypeScript/JavaScript), asks to add license headers to source files, needs
  to create or update a THIRD_PARTY_LICENSES file, or wants to ensure proper
  license compliance for derived works. Also trigger when the user provides
  a link to or local path of an upstream library and a folder of ported code
  and asks about licensing. Keywords: license, attribution, port, third-party,
  copyright, THIRD_PARTY_LICENSES, license header, BSD, MIT, Apache, derived
  work.
---

# License Updater

A skill for maintaining proper license attribution when incorporating
third-party code — especially code ported from one language to another.

## When to use

- User provides a **source library URL or local path** and a **folder with ported files**
- User asks to "add license headers", "update THIRD_PARTY_LICENSES", or
  "set up licensing for ported code"
- User mentions porting/translating code from an external project

## Inputs

The user provides two things:

1. **Source library URL or local path** — a GitHub (or similar) link to the
   upstream project whose code was ported, or a local filesystem path to a
   cloned/downloaded copy of the library.
   Examples: `https://github.com/sdwfrost/liblsoda`, `src/solver-tools/cvode-7.5.0`
2. **Ported files path** — a path (relative to the repo root) to the folder
   containing the ported TypeScript (or other language) files.
   Example: `src/solver-tools/lsoda`

The user may optionally provide:
- The **repo root URL** or local path to their own project
- The **target project's license** (if not, discover it from the repo)

## Workflow

Follow these steps in order:

### Step 1: Discover upstream license

If the user provided a **URL**, use `web_fetch` to open the source library URL.
If the user provided a **local path**, read files directly from the filesystem.

Look for:

1. A `LICENSE` or `COPYING` file in the repository root
2. License badges or mentions in the README
3. SPDX identifiers in `package.json`, `setup.py`, `Cargo.toml`, etc.
4. Copyright headers inside source files (especially the main source file)

Identify:
- **License type** (MIT, BSD-2-Clause, BSD-3-Clause, Apache-2.0, etc.)
- **Copyright holders** — every `Copyright (c)` line you find
- **Provenance chain** — if the upstream code itself was derived from earlier
  work, trace the chain (e.g., Fortran → C → C refactoring → shared library)

If the license is unclear, tell the user and ask for clarification before
proceeding.

### Step 2: Check compatibility

Verify that the upstream license is compatible with the target project's
license. Read `references/compatibility.md` for the compatibility matrix.

If there is a compatibility issue, **stop and inform the user** with a clear
explanation of the conflict and possible remedies.

### Step 3: Discover target project structure

If the user provided a repo URL, fetch it to understand the existing layout.
Look for:
- Existing `LICENSE` file
- Existing `THIRD_PARTY_LICENSES` or `NOTICES` file
- Existing copyright headers in files

### Step 4: Generate artifacts

Produce exactly three outputs. Use the script at
`scripts/generate_license_files.py` to generate them:

```bash
python3 /path/to/skill/scripts/generate_license_files.py \
  --license-type "<SPDX identifier, e.g. MIT>" \
  --copyright-holders "<semicolon-separated list>" \
  --project-name "<upstream project name>" \
  --project-url "<upstream project URL>" \
  --ported-files-path "<relative path in target repo>" \
  --target-copyright "<target project copyright line>" \
  --provenance "<semicolon-separated provenance chain>" \
  --output-dir "/home/claude/license-output"
```

The script creates three files in the output directory:

1. **`THIRD_PARTY_LICENSES`** — standalone file for the repo root documenting
   upstream license and provenance
2. **`HEADER_TEMPLATE.txt`** — comment block to prepend to each ported file
3. **`LICENSE_GUIDE.md`** — human-readable instructions explaining what to do

If the script is unavailable or doesn't fit the situation, generate the files
manually following the templates in `references/templates.md`.

### Step 5: Update the README License section

If the project has a `README.md` (or `README.MD`), update its **License**
section to mention the third-party code and link to `THIRD_PARTY_LICENSES`.

Rules:
- Look for an existing `## License` heading in the README
- If one exists, **append** a line about the ported code beneath the existing
  license statement — do not remove anything already there
- If there is no License section, add one at the end of the README (before
  any final "Contributing" or similar section)
- Keep it concise: one sentence identifying the ported component, its upstream
  source, and a link to the `THIRD_PARTY_LICENSES` file
- Use the template from `references/templates.md` → "README License section"

### Step 6: Present results

1. Give a brief summary:
   - What license was found
   - Whether it's compatible
   - What files were created or modified (THIRD_PARTY_LICENSES, headers,
     README)

## Important rules

- **Never modify the target project's main LICENSE file** unless the user
  explicitly asks. The main LICENSE covers the project as a whole.
- **Always preserve the full original license text** in THIRD_PARTY_LICENSES.
  Don't paraphrase or abbreviate license texts.
- **Include ALL copyright holders** from the upstream chain, not just the
  most recent one.
- **If the upstream repo contains multiple licenses** (e.g., dual-licensed),
  note all of them and ask the user which applies.
- If the target project already has a THIRD_PARTY_LICENSES file, **append**
  a new section rather than overwriting.
- The header template should be concise (under 15 lines) — point to
  THIRD_PARTY_LICENSES for full details rather than embedding the entire
  license in every file.
