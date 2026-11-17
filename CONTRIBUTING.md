## Contributing

### Generative AI Policy

First-time contributors should refrain from using generative AI without prior discussion. Unsolicited PRs that appear to rely significantly on AI may be closed immediately without review. AI tools must not be used in communication (with exception of translation).

### Getting started

Contributions fall a fork-pull request open source work flow:

1. Fork this repository
2. Clone your fork and build a development installation (see below)
3. Create a new branch (`git checkout -b your_new_branch_name`)
4. Commit changes on your branch
5. Push your branch to your fork
6. Submit pull request (PR) to the main `pyVBRc` repository

After PR submission, a number of automated tests will run and your changes
will be reviewed. Often this may involve some iteration, but we'll help you out
if you're new to git!

### Development installation

After cloning your fork, you can install `pyVBRc` with

```shell
pip install -e .[dev]
```

and then install all the development requirements with

```shell
pip install -r requirements_dev.txt
```

### running the test suite

After completing your development installation, you can run the test suite
with

```shell
pytest -v
```
To also generate a code coverage report

```shell
pytest -v --cov=./ --cov-report=html:coverage/
```
After tests run, you can then open the `coverage/index.html` file in your
preferred browser and check out where test coverage is poor.

Pull requests are required to maintain the code coverage of the code base, so
you'll have to write new tests to cover any new code. Bug fixes may not require
new tests if coverage is unchanged.

### code style

#### style checks
Code style is enforced with `pre-commit`. When you submit a pull request, the
`pre-commit.ci` bot will test whether your code passes the style requirements. If
that test fails, you can have `pre-commit.ci` autofix your pull request by adding
a comment to the pull request that reads

```
pre-commit.ci autofix
```

and the bot will try to fix the style failures. After it runs, if you still
have work to do on the PR, you'll want to run

```
git pull
```

locally so that you can pull in the bot's changes to your branch.

#### using pre-commit during development

If you'd like to ensure that your code passes all style requirements while
you're developing, you can use pre-commit throughout development. After
building your development installation, you can initialize pre-commit within
your repository with

```shell
pre-commit install
```

and your code changes will be checked every time you commit.
