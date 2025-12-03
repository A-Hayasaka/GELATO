# GitHub Pages Setup for GELATO Documentation

## Setup Instructions

### 1. Enable GitHub Pages in Repository Settings

1. Go to your GitHub repository: `https://github.com/istellartech/GELATO`
2. Click on **Settings** tab
3. In the left sidebar, click **Pages**
4. Under **Build and deployment**:
   - **Source**: Select "GitHub Actions"
5. Save the settings

### 2. Push the Changes

```bash
git add .github/workflows/deploy-docs.yml docs/.nojekyll
git commit -m "Add GitHub Pages deployment workflow"
git push origin add_document
```

### 3. Merge to Main Branch

Once the PR is merged to `master` or `main` branch, the documentation will automatically build and deploy.

### 4. Access the Documentation

After deployment completes (usually 2-5 minutes), the documentation will be available at:

```
https://istellartech.github.io/GELATO/
```

## Automatic Updates

The documentation will automatically rebuild and redeploy whenever:
- Changes are pushed to the `master` or `main` branch
- Manual workflow trigger via GitHub Actions tab

## Local Testing

To test the documentation build locally:

```bash
cd docs
make html
```

Then open `docs/_build/html/index.html` in your browser.

## Badges for README

Add these badges to your README.md:

```markdown
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://istellartech.github.io/GELATO/)
[![Build Status](https://github.com/istellartech/GELATO/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/istellartech/GELATO/actions/workflows/deploy-docs.yml)
```

## Troubleshooting

### Build Fails

1. Check the Actions tab for error messages
2. Ensure all dependencies are in `docs/requirements.txt`
3. Test locally with `cd docs && make html`

### Missing C++ Documentation

- Ensure Doxygen is installed (handled by workflow)
- Check that `Doxyfile` is properly configured in `docs/`

### Broken Links or Missing Assets

- Ensure `.nojekyll` file is present
- Check that all paths in Sphinx configuration are relative
- Verify `html_static_path` in `conf.py`
