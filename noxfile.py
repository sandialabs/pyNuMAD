import nox

@nox.session
def tests(session):
    session.install(".")
    session.install("pytest")
    session.run("pytest")

@nox.session
def lint(session):
    session.install("flake8")
    session.run("flake8", "--import-order-style", "google")
    
@nox.session
def docs(session):
    session.install(".")
    session.install("sphinx")
    session.install("sphinx_rtd_theme")
    session.install("sphinxcontrib-bibtex")
    session.run("sphinx-build", "-M", "html", "docs/", "docs/_build")

@nox.session
def serve(session):
    session.run("python", "-m", "http.server", "-b", "localhost", "-d", "docs/_build/html", "8085")