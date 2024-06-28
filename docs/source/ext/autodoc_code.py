import inspect

def setup(app):
    app.connect('autodoc-process-docstring', append_source_code)
    return {'version': '0.1'}

def append_source_code(app, what, name, obj, options, lines):
    if what == "function":
        try:
            source = inspect.getsource(obj)
            lines.append(".. code-block:: python\n")
            lines.extend(f"    {line}" for line in source.splitlines())
        except Exception as e:
            lines.append(".. code-block:: python\n")
            lines.append(f"    # Error retrieving source code: {str(e)}")
