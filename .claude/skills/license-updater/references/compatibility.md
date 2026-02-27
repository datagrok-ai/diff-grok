# License Compatibility Matrix

Quick reference for whether an upstream license is compatible with common
target project licenses. "Compatible" means the ported code can be included
in a project under the target license.

## Permissive → Permissive (almost always OK)

| Upstream        | Target: MIT | Target: BSD-2 | Target: BSD-3 | Target: Apache-2.0 |
|-----------------|:-----------:|:-------------:|:-------------:|:------------------:|
| MIT             | ✅          | ✅            | ✅            | ✅                 |
| BSD-2-Clause    | ✅          | ✅            | ✅            | ✅                 |
| BSD-3-Clause    | ✅          | ✅            | ✅            | ✅                 |
| Apache-2.0      | ✅*         | ✅*           | ✅*           | ✅                 |
| ISC             | ✅          | ✅            | ✅            | ✅                 |
| Zlib            | ✅          | ✅            | ✅            | ✅                 |
| Unlicense       | ✅          | ✅            | ✅            | ✅                 |

*Apache-2.0 has a patent grant clause. When including Apache-2.0 code in a
non-Apache project, the patent grant still applies to that code. This is
generally fine but worth noting.

## Copyleft → Permissive (usually NOT OK)

| Upstream        | Target: MIT | Target: BSD | Target: Apache-2.0 |
|-----------------|:-----------:|:-----------:|:------------------:|
| GPL-2.0-only    | ❌          | ❌          | ❌                 |
| GPL-3.0-only    | ❌          | ❌          | ❌                 |
| LGPL-2.1        | ⚠️*         | ⚠️*         | ⚠️*                |
| LGPL-3.0        | ⚠️*         | ⚠️*         | ⚠️*                |
| MPL-2.0         | ⚠️**        | ⚠️**        | ⚠️**               |
| AGPL-3.0        | ❌          | ❌          | ❌                 |

*LGPL allows linking but not incorporating source code into a permissive
project without the LGPL applying to the modified files.

**MPL-2.0 is file-level copyleft: the ported files themselves must stay
under MPL-2.0, but the rest of the project can be permissive.

## What "compatible" means in practice

For permissive licenses (MIT, BSD, Apache, ISC, Zlib):
- Keep the original copyright notice and license text
- Include them in a THIRD_PARTY_LICENSES file
- Add attribution headers to ported files
- The overall project can use any license

For copyleft (GPL, LGPL, AGPL):
- Generally cannot include in a permissive project
- The copyleft license would need to apply to the derivative work
- STOP and inform the user about the conflict

## Special conditions by license type

### BSD-3-Clause
Extra requirement: cannot use the names of copyright holders or contributors
to endorse or promote the derived product without written permission.

### Apache-2.0
- Includes explicit patent grant
- Requires NOTICE file preservation if one exists upstream
- Must state significant changes made to ported files

### MIT
Minimal requirements: preserve copyright notice and license text.

### ISC
Functionally equivalent to MIT. Same requirements.
