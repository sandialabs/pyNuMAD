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
    session.cd("docs/")
    session.run("make html")

@nox.session
def serve(session):
    session.run("python", "-m", "http.server", "-b", "localhost", "-d", "docs/_build/html", "8085")
   
@nox.session
def check_style(session):
    session.install("black")
    session.run("black", "--check")
    
@nox.session
def enforce_style(session):
    session.install("black")
    session.run("black", "src")