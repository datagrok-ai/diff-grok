#!/usr/bin/env python3
"""Generate license attribution files for ported third-party code.

This script produces three files:
  1. THIRD_PARTY_LICENSES — for the repo root
  2. HEADER_TEMPLATE.txt  — comment block to prepend to each ported file
  3. LICENSE_GUIDE.md      — human-readable instructions

Usage:
    python3 generate_license_files.py \
        --license-type "MIT" \
        --copyright-holders "Copyright (c) 2011 Foo;Copyright (c) 2020 Bar" \
        --project-name "liblsoda" \
        --project-url "https://github.com/sdwfrost/liblsoda" \
        --ported-files-path "src/solver-tools/lsoda" \
        --target-copyright "Copyright (c) 2025 Datagrok" \
        --output-dir "./license-output"

Optional:
    --provenance       Semicolon-separated provenance chain
    --description      Brief description of what the ported code does
    --license-text     Full license text (if not a well-known SPDX type)
    --target-license   Target project license type (default: MIT)
    --file-ext         File extension for header comment style (default: ts)
    --append           If set, append to existing THIRD_PARTY_LICENSES
"""

import argparse
import os
import sys
import textwrap

# Well-known license texts keyed by SPDX identifier
KNOWN_LICENSES = {
    "MIT": textwrap.dedent("""\
        The MIT License

        {COPYRIGHT}

        Permission is hereby granted, free of charge, to any person obtaining
        a copy of this software and associated documentation files (the
        "Software"), to deal in the Software without restriction, including
        without limitation the rights to use, copy, modify, merge, publish,
        distribute, sublicense, and/or sell copies of the Software, and to
        permit persons to whom the Software is furnished to do so, subject to
        the following conditions:

        The above copyright notice and this permission notice shall be
        included in all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
        EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
        MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
        NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
        LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
        OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
        WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."""),

    "BSD-2-Clause": textwrap.dedent("""\
        BSD 2-Clause License

        {COPYRIGHT}

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

        1. Redistributions of source code must retain the above copyright notice,
           this list of conditions and the following disclaimer.

        2. Redistributions in binary form must reproduce the above copyright notice,
           this list of conditions and the following disclaimer in the documentation
           and/or other materials provided with the distribution.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
        AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
        IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
        FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
        DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
        SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
        CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
        OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
        OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."""),

    "BSD-3-Clause": textwrap.dedent("""\
        BSD 3-Clause License

        {COPYRIGHT}

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

        1. Redistributions of source code must retain the above copyright notice,
           this list of conditions and the following disclaimer.

        2. Redistributions in binary form must reproduce the above copyright notice,
           this list of conditions and the following disclaimer in the documentation
           and/or other materials provided with the distribution.

        3. Neither the name of the copyright holder nor the names of its
           contributors may be used to endorse or promote products derived from
           this software without specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
        AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
        IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
        FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
        DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
        SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
        CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
        OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
        OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."""),

    "Apache-2.0": textwrap.dedent("""\
        Apache License, Version 2.0

        {COPYRIGHT}

        Licensed under the Apache License, Version 2.0 (the "License");
        you may not use this file except in compliance with the License.
        You may obtain a copy of the License at

            http://www.apache.org/licenses/LICENSE-2.0

        Unless required by applicable law or agreed to in writing, software
        distributed under the License is distributed on an "AS IS" BASIS,
        WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
        See the License for the specific language governing permissions and
        limitations under the License."""),

    "ISC": textwrap.dedent("""\
        ISC License

        {COPYRIGHT}

        Permission to use, copy, modify, and/or distribute this software for any
        purpose with or without fee is hereby granted, provided that the above
        copyright notice and this permission notice appear in all copies.

        THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
        WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
        ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
        WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
        ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
        OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE."""),
}

# Compatibility matrix: (upstream, target) -> compatible?
# True = compatible, False = incompatible, None = conditional
COMPAT = {
    ("MIT", "MIT"): True,
    ("MIT", "BSD-2-Clause"): True,
    ("MIT", "BSD-3-Clause"): True,
    ("MIT", "Apache-2.0"): True,
    ("BSD-2-Clause", "MIT"): True,
    ("BSD-2-Clause", "BSD-2-Clause"): True,
    ("BSD-2-Clause", "BSD-3-Clause"): True,
    ("BSD-2-Clause", "Apache-2.0"): True,
    ("BSD-3-Clause", "MIT"): True,
    ("BSD-3-Clause", "BSD-2-Clause"): True,
    ("BSD-3-Clause", "BSD-3-Clause"): True,
    ("BSD-3-Clause", "Apache-2.0"): True,
    ("Apache-2.0", "MIT"): True,
    ("Apache-2.0", "BSD-2-Clause"): True,
    ("Apache-2.0", "BSD-3-Clause"): True,
    ("Apache-2.0", "Apache-2.0"): True,
    ("ISC", "MIT"): True,
    ("ISC", "BSD-2-Clause"): True,
    ("ISC", "BSD-3-Clause"): True,
    ("ISC", "Apache-2.0"): True,
}

COMMENT_STYLES = {
    "ts": ("/**", " *", " */"),
    "js": ("/**", " *", " */"),
    "tsx": ("/**", " *", " */"),
    "jsx": ("/**", " *", " */"),
    "py": ("#", "#", "#"),
    "c": ("/*", " *", " */"),
    "cpp": ("/*", " *", " */"),
    "h": ("/*", " *", " */"),
    "rs": ("//!", "//!", "//!"),
    "go": ("//", "//", "//"),
}


def parse_list(value: str) -> list[str]:
    """Split a semicolon-separated string into a list of stripped items."""
    if not value:
        return []
    return [item.strip() for item in value.split(";") if item.strip()]


def get_license_text(license_type: str, copyright_holders: list[str],
                     custom_text: str | None = None) -> str:
    """Return the full license text with copyright holders inserted."""
    if custom_text:
        return custom_text

    template = KNOWN_LICENSES.get(license_type)
    if not template:
        return (f"[License text for {license_type} not available.\n"
                f"Please paste the full license text here manually.]")

    copyright_block = "\n".join(copyright_holders)
    return template.replace("{COPYRIGHT}", copyright_block)


def check_compatibility(upstream: str, target: str) -> tuple[bool, str]:
    """Check license compatibility and return (ok, message)."""
    key = (upstream, target)
    if key in COMPAT:
        if COMPAT[key]:
            return True, f"{upstream} is compatible with {target}."
        else:
            return False, (f"⚠️  {upstream} is NOT compatible with {target}. "
                           f"The ported code cannot be included under {target}.")

    # Assume copyleft is incompatible with permissive
    copyleft = {"GPL-2.0-only", "GPL-3.0-only", "AGPL-3.0-only",
                "GPL-2.0-or-later", "GPL-3.0-or-later", "AGPL-3.0-or-later"}
    if upstream in copyleft:
        return False, (f"⚠️  {upstream} is copyleft and cannot be included "
                       f"in a {target} project without the entire project "
                       f"adopting {upstream}.")

    return True, (f"Compatibility between {upstream} and {target} is not in "
                  f"the built-in matrix. Manual review recommended.")


def generate_third_party(args, copyright_holders: list[str],
                         provenance: list[str]) -> str:
    """Generate the THIRD_PARTY_LICENSES content."""
    lines = []
    lines.append("Third-Party Licenses")
    lines.append("=" * 20)
    lines.append("")
    lines.append("This file documents the licenses for third-party code included in or adapted")
    lines.append("by this project.")
    lines.append("")
    lines.append("")
    lines.append("=" * 80)
    lines.append(f"{args.project_name}")
    lines.append("=" * 80)
    lines.append("")
    lines.append(f"Source:       {args.project_url}")
    lines.append(f"License:      {args.license_type}")
    lines.append(f"Files:        {args.ported_files_path}")
    lines.append("")

    if args.description:
        lines.append(args.description)
        lines.append("")

    if provenance:
        lines.append("Provenance chain:")
        lines.append("")
        for i, step in enumerate(provenance, 1):
            lines.append(f"  {i}. {step}")
        lines.append("")

    license_text = get_license_text(
        args.license_type, copyright_holders, args.license_text
    )
    lines.append(license_text)
    lines.append("")

    return "\n".join(lines)


def generate_header(args, copyright_holders: list[str]) -> str:
    """Generate the header template for ported files."""
    ext = args.file_ext
    start, mid, end = COMMENT_STYLES.get(ext, ("/**", " *", " */"))

    lines = []
    lines.append(start)

    if args.description:
        lines.append(f"{mid} {args.description}")
        lines.append(f"{mid}")

    lines.append(f"{mid} Port of {args.project_name} ({args.project_url})")
    lines.append(f"{mid}")

    for holder in copyright_holders:
        lines.append(f"{mid} {holder}")
    if args.target_copyright:
        lines.append(f"{mid} {args.target_copyright}")

    lines.append(f"{mid}")
    lines.append(f"{mid} Licensed under the {args.license_type} License.")
    lines.append(f"{mid} See THIRD_PARTY_LICENSES for the full license text and provenance.")
    lines.append(end)
    lines.append("")

    return "\n".join(lines)


def generate_guide(args, copyright_holders: list[str],
                   compat_ok: bool, compat_msg: str) -> str:
    """Generate the LICENSE_GUIDE.md content."""
    status = "✅ Compatible" if compat_ok else "❌ Incompatible"
    target_lic = args.target_license

    lines = []
    lines.append(f"# License Attribution Guide for {args.project_name}")
    lines.append("")
    lines.append("## File structure")
    lines.append("")
    lines.append("```")
    lines.append("your-repo/")
    lines.append("├── LICENSE                      ← do not modify")
    lines.append("├── THIRD_PARTY_LICENSES         ← add/update this file")

    # Show the ported files path
    parts = args.ported_files_path.split("/")
    indent = "│   "
    for i, part in enumerate(parts[:-1]):
        lines.append(f"{indent}└── {part}/")
        indent += "    "
    lines.append(f"{indent}└── {parts[-1]}/")
    lines.append(f"{indent}    └── *.{args.file_ext}          ← add header to each file")

    lines.append("```")
    lines.append("")
    lines.append("## Steps")
    lines.append("")
    lines.append("### 1. Add `THIRD_PARTY_LICENSES` to the repository root")
    lines.append("")

    if args.append:
        lines.append("Append the provided section to your existing `THIRD_PARTY_LICENSES` file.")
    else:
        lines.append("Place the provided `THIRD_PARTY_LICENSES` file next to your existing `LICENSE`.")

    lines.append("")
    lines.append(f"### 2. Add header to each `.{args.file_ext}` file in `{args.ported_files_path}/`")
    lines.append("")
    lines.append("Prepend the content of `HEADER_TEMPLATE.txt` to the top of every")
    lines.append(f"`.{args.file_ext}` file in the `{args.ported_files_path}/` directory.")
    lines.append("")
    lines.append("### 3. `LICENSE` — no changes needed")
    lines.append("")
    lines.append(f"Your project license ({target_lic}) stays as-is.")
    lines.append("`THIRD_PARTY_LICENSES` handles attribution for ported code separately.")
    lines.append("")
    lines.append("### 4. README — optional")
    lines.append("")
    lines.append("Consider adding a note in your README's License section:")
    lines.append("")
    lines.append("```markdown")
    lines.append("## License")
    lines.append("")
    lines.append(f"This project is licensed under the [{target_lic} License](LICENSE).")
    lines.append("")
    lines.append(f"Some components are ported from [{args.project_name}]({args.project_url}).")
    lines.append("See [THIRD_PARTY_LICENSES](THIRD_PARTY_LICENSES) for details.")
    lines.append("```")
    lines.append("")
    lines.append("## Compatibility check")
    lines.append("")
    lines.append(f"- Upstream license: **{args.license_type}**")
    lines.append(f"- Target license:   **{target_lic}**")
    lines.append(f"- Status:           {status}")
    lines.append("")
    lines.append(compat_msg)

    # Special notes
    if args.license_type == "BSD-3-Clause":
        lines.append("")
        lines.append("### BSD-3-Clause special condition")
        lines.append("")
        lines.append("You may NOT use the names of the original copyright holders")
        lines.append("to endorse or promote your product without their written permission.")

    if args.license_type == "Apache-2.0":
        lines.append("")
        lines.append("### Apache-2.0 special conditions")
        lines.append("")
        lines.append("- The patent grant from Apache-2.0 applies to the ported code")
        lines.append("- If the upstream repo has a NOTICE file, its contents must be preserved")
        lines.append("- You must state significant changes made to the ported files")

    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Generate license attribution files for ported code."
    )
    parser.add_argument("--license-type", required=True,
                        help="SPDX license identifier (e.g., MIT, BSD-3-Clause)")
    parser.add_argument("--copyright-holders", required=True,
                        help="Semicolon-separated copyright lines")
    parser.add_argument("--project-name", required=True,
                        help="Name of the upstream project")
    parser.add_argument("--project-url", required=True,
                        help="URL to the upstream project")
    parser.add_argument("--ported-files-path", required=True,
                        help="Relative path in target repo to ported files")
    parser.add_argument("--target-copyright", default="",
                        help="Copyright line for the target project")
    parser.add_argument("--target-license", default="MIT",
                        help="Target project license (default: MIT)")
    parser.add_argument("--provenance", default="",
                        help="Semicolon-separated provenance chain")
    parser.add_argument("--description", default="",
                        help="Brief description of the ported code")
    parser.add_argument("--license-text", default=None,
                        help="Full custom license text (if not a known SPDX type)")
    parser.add_argument("--file-ext", default="ts",
                        help="File extension for comment style (default: ts)")
    parser.add_argument("--output-dir", default="./license-output",
                        help="Output directory for generated files")
    parser.add_argument("--append", action="store_true",
                        help="Append to existing THIRD_PARTY_LICENSES")

    args = parser.parse_args()

    copyright_holders = parse_list(args.copyright_holders)
    provenance = parse_list(args.provenance)

    if not copyright_holders:
        print("Error: at least one copyright holder is required.", file=sys.stderr)
        sys.exit(1)

    # Check compatibility
    compat_ok, compat_msg = check_compatibility(args.license_type, args.target_license)
    if not compat_ok:
        print(f"\n⚠️  LICENSE COMPATIBILITY WARNING:\n{compat_msg}\n", file=sys.stderr)
        print("Files will still be generated, but review carefully before using.",
              file=sys.stderr)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate files
    third_party = generate_third_party(args, copyright_holders, provenance)
    header = generate_header(args, copyright_holders)
    guide = generate_guide(args, copyright_holders, compat_ok, compat_msg)

    # Write files
    tp_path = os.path.join(args.output_dir, "THIRD_PARTY_LICENSES")
    with open(tp_path, "w") as f:
        f.write(third_party)
    print(f"✅ Created: {tp_path}")

    header_path = os.path.join(args.output_dir, "HEADER_TEMPLATE.txt")
    with open(header_path, "w") as f:
        f.write(header)
    print(f"✅ Created: {header_path}")

    guide_path = os.path.join(args.output_dir, "LICENSE_GUIDE.md")
    with open(guide_path, "w") as f:
        f.write(guide)
    print(f"✅ Created: {guide_path}")

    print(f"\n{'✅' if compat_ok else '⚠️'}  License compatibility: {compat_msg}")
    print(f"\nAll files written to: {args.output_dir}/")


if __name__ == "__main__":
    main()
