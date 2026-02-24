import nox

@nox.session
def tests(session):
    """Run tests."""
    session.install(".")
    session.install("pytest")
    session.run("pytest")

@nox.session
def lint(session):
    """Lint."""
    session.install("flake8")
    session.run("flake8", "--import-order-style", "google")
    
@nox.session
def docs(session):
    """Generate documentation."""
    session.run("pip", "install", "-e", ".")
    session.install("sphinx")
    session.install("pydata-sphinx-theme")
    session.install("sphinxcontrib-bibtex")
    session.cd("docs/")
    session.run("make", "html")

@nox.session
def serve(session):
    """Serve documentation. Port can be specified as a positional argument."""
    try:
        port = session.posargs[0]
    except IndexError:
        port = "8085"
    session.run("python", "-m", "http.server", "-b", "localhost", "-d", "docs/_build/html", port)
   
@nox.session
def check_style(session):
    """Check if code follows black style."""
    session.install("black")
    session.run("black", "--check", "src")
    
@nox.session
def enforce_style(session):
    """Apply black style to code base."""
    session.install("black")
    session.run("black", "src")
