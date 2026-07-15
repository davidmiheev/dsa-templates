# Rust + Sphinx Integration Proposal

## Current Problems

1. **Rust code is duplicated** — manually copy-pasted into `.rst` files as `.. code-block:: rust`, easily goes out of sync
2. **Rust modules mix impl + `fn main()`** — files include test scaffolding that shouldn't appear in docs
3. **No unified view** — docs are organized by language (Python section + Rust section) rather than by topic
4. **`docs/source/_static/` is not in git** — logo was lost on rebuild, only fixed after manually staging

---

## Proposed Improvements

### 1. Use `literalinclude` Instead of Hardcoded Code Blocks

Replace manually maintained code in `.rst` files with `.. literalinclude::` directives that read directly from source:

```rst
Rust Implementation
-------------------

.. literalinclude:: ../../sorted_list/sorted_list.rs
   :language: rust
   :linenos:
   :dedent: 4
   :start-after: pub struct SortedList
```

**Pros**: docs stay in sync with source automatically
**Cons**: if source has test scaffolding, it bleeds into docs

### 2. Separate Library from Examples

Split each Rust module into a clean library + test driver:

```
sorted_list/
├── src/
│   ├── lib.rs          # clean impl only (for docs)
│   └── main.rs         # test/example driver
└── sorted_list.rs      # (remove, moved to src/)
```

Same for Python — create a parallel `src/` structure or keep clean impls in the main `.py` file.

### 3. Recommended Final Structure

```
sorted_list/
├── __init__.py         # from .sorted_list import SortedList
├── sorted_list.py      # clean Python impl
├── src/
│   ├── lib.rs          # clean Rust impl
│   └── main.rs         # test driver (not in docs)
```

- Python: use existing `.. automodule::` for docstrings
- Rust: use `.. literalinclude::` from `src/lib.rs`

### 4. Fix Git Static Files

Ensure `docs/source/_static/logo.png` is committed so ReadTheDocs can use it.

---

## Implementation Roadmap

### Phase 1: Quick Win (low effort, high value)
- [ ] Convert all `.rst` files to use `.. literalinclude::` for Rust code
- [ ] Commit `docs/source/_static/` to git
- [ ] Test that docs rebuild correctly from source

### Phase 2: Restructure (medium effort)
- [ ] Create `src/lib.rs` + `src/main.rs` in each Rust topic
- [ ] Clean up Python files (move test scaffolding out)
- [ ] Update all `literalinclude` paths

### Phase 3: Polish (optional)
- [ ] Add tabbed code blocks (Python / Rust) per topic using `sphinx-tabs`
- [ ] Add custom Sphinx extension to parse Rust `///` doc comments
- [ ] Generate API reference for Rust from `lib.rs`

---

## Commands to Know

```bash
# Build docs
cd docs/source && uv run sphinx-build -b html . _build/html

# Serve locally
cd docs/source && uv run python -m http.server 8000 -d _build/html

# ReadTheDocs config
cat .readthedocs.yaml
```
