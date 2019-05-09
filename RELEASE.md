
Make sure packaging tools are up to date.

```
. venv/bin/activate
pip install --upgrade pip
pip install --upgrade twine
pip install --upgrade setuptools wheel
```

Bump version; document changes.

```
git tag v0.1.7
vim cathpy/__init__.py
vim Changelog
```

Package release.

```
python3 setup.py sdist bdist_wheel
```

Upload files.

```
python3 -m twine upload dist/*
```