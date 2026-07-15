"""Sphinx extension that integrates ``cargo doc`` output into the build.

Adds a reST directive ``.. rust-api:: <crate>`` that renders as a plain
heading + paragraph + link, identical in style to any other RTD section.
The directive's link points to ``_rust_doc/<crate>/index.html``.

The entire ``target/doc/`` tree is copied into the Sphinx output dir under
``_rust_doc/`` via the ``build-finished`` hook, so all the cargo doc
inter-linking (crates sidebar, search, fonts) keeps working.

If ``cargo`` is not on PATH at build time, the directive still renders a
link; a warning is logged but the build does not fail.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

from docutils import nodes
from docutils.parsers.rst import Directive
from sphinx.application import Sphinx
from sphinx.util import logging

logger = logging.getLogger(__name__)


# Repository layout assumptions (resolved from the Sphinx project's root):
#   conf.py is at <repo>/docs/source/conf.py
#   workspace Cargo.toml is at <repo>/Cargo.toml
#   cargo doc output is at <repo>/target/doc/
REPO_ROOT = Path(__file__).resolve().parents[3]
CARGO_TARGET_DOC = REPO_ROOT / "target" / "doc"


def _cargo_available() -> bool:
    return shutil.which("cargo") is not None


def workspace_member_dirs() -> list[str]:
    """Read ``members`` from the workspace ``Cargo.toml`` at REPO_ROOT.

    Used as a cheap "is cargo doc output fresh?" probe — if any member
    crate's ``index.html`` already exists under ``target/doc/``, we skip
    rebuilding. Returns an empty list on any parse failure (the build
    then falls through to running ``cargo doc``).
    """
    cargo_toml = REPO_ROOT / "Cargo.toml"
    if not cargo_toml.exists():
        return []
    try:
        import tomllib  # py3.11+
    except ImportError:
        try:
            import tomli as tomllib  # type: ignore
        except ImportError:
            return []
    try:
        data = tomllib.loads(cargo_toml.read_text())
        return list(data.get("workspace", {}).get("members", []))
    except Exception:
        return []


def _ensure_cargo_doc_built(app: Sphinx) -> None:
    """Run ``cargo doc --no-deps --workspace`` once at build start if needed.

    Skipped if cargo is missing or output is already current.
    """
    if not _cargo_available():
        logger.warning(
            "rustdoc: cargo not found on PATH; skipping `cargo doc`. "
            "If `_rust_doc/` is empty, install Rust or run "
            "`cargo doc --no-deps --workspace` manually before building docs."
        )
        return

    cargo_toml = REPO_ROOT / "Cargo.toml"
    if not cargo_toml.exists():
        logger.warning("rustdoc: %s not found; nothing to document.", cargo_toml)
        return

    docs_dir = CARGO_TARGET_DOC
    docs_dir_exists = docs_dir.is_dir()

    # We don't emit a top-level target/doc/index.html (no workspace src/lib.rs),
    # so use any crate's index.html as a "did cargo doc run?" proxy.
    if docs_dir_exists and any((docs_dir / c / "index.html").exists()
                               for c in workspace_member_dirs()):
        logger.info("rustdoc: `cargo doc` already built at %s.", docs_dir)
        return

    logger.info("rustdoc: running `cargo doc --no-deps --workspace` ... this may take a while.")
    try:
        subprocess.run(
            ["cargo", "doc", "--no-deps", "--workspace", "--document-private-items"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=600,
        )
    except subprocess.CalledProcessError as e:
        logger.warning(
            "rustdoc: `cargo doc` failed (exit %s). Rust API pages will be broken.\n%s",
            e.returncode,
            (e.stderr or e.stdout or "")[:2000],
        )
    except subprocess.TimeoutExpired:
        logger.warning("rustdoc: `cargo doc` timed out after 600s.")


def _copy_cargo_doc_into_outdir(app: Sphinx, exception: Exception | None) -> None:
    """Copy ``target/doc/`` tree into ``<outdir>/_rust_doc/``.

    Fires on the ``build-finished`` hook, which runs *after* all pages are
    written and *after* Sphinx is done mutating ``outdir``. The
    ``exception`` argument is non-None if the build itself failed; we skip
    copying in that case so we don't ship stale docs alongside broken ones.
    """
    if exception is not None:
        return

    outdir = Path(app.outdir)
    dest = outdir / "_rust_doc"

    # Clean up any previous run so we don't leave stale files behind.
    if dest.exists():
        shutil.rmtree(dest)

    if not CARGO_TARGET_DOC.is_dir():
        logger.warning("rustdoc: %s not found; skipping copy.", CARGO_TARGET_DOC)
        return

    # Copy the *entire* cargo doc tree as one unit so cross-crate links
    # (sidebar, search index, crates.js, static files) keep working.
    shutil.copytree(CARGO_TARGET_DOC, dest)
    logger.info("rustdoc: copied %s -> %s", CARGO_TARGET_DOC, dest)


def _make_reference(text: str, refuri: str) -> nodes.reference:
    """Build a docutils reference node that actually renders as <a href>.

    The trick is that docutils infers the ``rawsource`` from the text
    children when it serialises. Passing ``rawsource=""`` forces the whole
    node to be treated as escaped text, which is what produced the
    ``&lt;reference ...&gt;`` we were seeing. Letting docutils compute it
    from the text child fixes the rendering.
    """
    node = nodes.reference(refuri=refuri, internal=False)
    node += nodes.Text(text, text)
    return node


class RustApiDirective(Directive):
    """Render an "API Documentation" section that links to a crate's cargo doc.

    Usage::

        .. rust-api:: fenwick
           :title: Fenwick Tree (optional; defaults to crate name)
           :description: Short one-line description (optional)

    The output is a regular ``<h2>`` heading + paragraph + link, identical
    in style to prose anywhere else on the page.
    """

    required_arguments = 1  # the crate name
    optional_arguments = 0
    has_content = False
    option_spec = {
        "title": lambda s: s.strip(),
        "description": lambda s: s.strip(),
    }

    def run(self) -> list[nodes.Node]:
        crate = self.arguments[0]
        title = self.options.get("title", crate)
        description = self.options.get(
            "description",
            f"Generated API reference for the {crate} crate.",
        )

        # Return a flat list of nodes (heading, description paragraph, link
        # paragraph) so docutils drops them into the surrounding document
        # at the directive's location. Returning a nested ``section`` would
        # make this a sub-section under whatever <section> the page is in,
        # which is what we want for headings that nest under the parent
        # heading level. We use raw heading + paragraphs to keep styling
        # identical to other prose on the page.
        result: list[nodes.Node] = []

        heading = nodes.title(
            "",
            f"{title} \u2014 Rust API",
        )
        result.append(heading)

        if description:
            result.append(nodes.paragraph("", description))

        link_para = nodes.paragraph(
            "",
            "",
            _make_reference(
                "View the rendered API reference \u2192",
                _rust_doc_url(crate),
            ),
        )
        result.append(link_para)

        return result


def _rust_doc_url(crate: str) -> str:
    """Resolve a cargo-doc URL relative to a page under Sphinx's ``outdir``.

    Sphinx pages live at ``<outdir>/<page>.html``; cargo doc lives at
    ``<outdir>/_rust_doc/<crate>/``. From any topic page, the link is
    ``_rust_doc/<crate>/index.html`` (relative).
    """
    return f"_rust_doc/{crate}/index.html"


def setup(app: Sphinx) -> dict[str, Any]:
    app.add_directive("rust-api", RustApiDirective)
    app.connect("builder-inited", _ensure_cargo_doc_built)
    app.connect("build-finished", _copy_cargo_doc_into_outdir)
    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
