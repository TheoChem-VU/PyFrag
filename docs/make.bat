:: filepath: c:\Users\Siebe\OneDrive - Vrije Universiteit Amsterdam\PhD\Scripting\local_packages\PyFrag\docs\make.bat
@echo off
REM Windows batch file for building Sphinx documentation

set SPHINXBUILD=sphinx-build
set BUILDDIR=build

if "%1" == "" (
    echo Please specify a target. Available targets are:
    echo   html       to make standalone HTML files
    echo   dirhtml    to make HTML files named index.html in directories
    echo   singlehtml to make a single large HTML file
    echo   pickle     to make pickle files
    echo   json       to make JSON files
    echo   htmlhelp   to make HTML files and a HTML help project
    echo   qthelp     to make HTML files and a qthelp project
    echo   epub       to make an epub
    echo   latex      to make LaTeX files
    echo   latexpdf   to make LaTeX files and run them through pdflatex
    echo   text       to make text files
    echo   man        to make manual pages
    echo   gettext    to make PO message catalogs
    echo   changes    to make an overview of all changed/added/deprecated items
    echo   linkcheck  to check all external links for integrity
    echo   doctest    to run all doctests embedded in the documentation
    echo   clean      to clean up the build directory
    exit /b 1
)

if "%1" == "clean" (
    echo Cleaning up build directory...
    rmdir /s /q %BUILDDIR%
    exit /b 0
)

if "%1" == "html" (
    %SPHINXBUILD% -b html . %BUILDDIR%/html
    echo Build finished. The HTML pages are in %BUILDDIR%/html.
    exit /b 0
)

if "%1" == "dirhtml" (
    %SPHINXBUILD% -b dirhtml . %BUILDDIR%/dirhtml
    echo Build finished. The HTML pages are in %BUILDDIR%/dirhtml.
    exit /b 0
)

if "%1" == "singlehtml" (
    %SPHINXBUILD% -b singlehtml . %BUILDDIR%/singlehtml
    echo Build finished. The HTML page is in %BUILDDIR%/singlehtml.
    exit /b 0
)

if "%1" == "pickle" (
    %SPHINXBUILD% -b pickle . %BUILDDIR%/pickle
    echo Build finished. The pickle files are in %BUILDDIR%/pickle.
    exit /b 0
)

if "%1" == "json" (
    %SPHINXBUILD% -b json . %BUILDDIR%/json
    echo Build finished. The JSON files are in %BUILDDIR%/json.
    exit /b 0
)

if "%1" == "htmlhelp" (
    %SPHINXBUILD% -b htmlhelp . %BUILDDIR%/htmlhelp
    echo Build finished. The HTML help project is in %BUILDDIR%/htmlhelp.
    exit /b 0
)

if "%1" == "qthelp" (
    %SPHINXBUILD% -b qthelp . %BUILDDIR%/qthelp
    echo Build finished. The Qt help project is in %BUILDDIR%/qthelp.
    exit /b 0
)

if "%1" == "epub" (
    %SPHINXBUILD% -b epub . %BUILDDIR%/epub
    echo Build finished. The epub file is in %BUILDDIR%/epub.
    exit /b 0
)

if "%1" == "latex" (
    %SPHINXBUILD% -b latex . %BUILDDIR%/latex
    echo Build finished. The LaTeX files are in %BUILDDIR%/latex.
    exit /b 0
)

if "%1" == "latexpdf" (
    %SPHINXBUILD% -b latex . %BUILDDIR%/latex
    echo Running LaTeX files through pdflatex...
    pushd %BUILDDIR%/latex
    pdflatex *.tex
    popd
    echo pdflatex finished. The PDF files are in %BUILDDIR%/latex.
    exit /b 0
)

if "%1" == "text" (
    %SPHINXBUILD% -b text . %BUILDDIR%/text
    echo Build finished. The text files are in %BUILDDIR%/text.
    exit /b 0
)

if "%1" == "man" (
    %SPHINXBUILD% -b man . %BUILDDIR%/man
    echo Build finished. The manual pages are in %BUILDDIR%/man.
    exit /b 0
)

if "%1" == "gettext" (
    %SPHINXBUILD% -b gettext . %BUILDDIR%/locale
    echo Build finished. The message catalogs are in %BUILDDIR%/locale.
    exit /b 0
)

if "%1" == "changes" (
    %SPHINXBUILD% -b changes . %BUILDDIR%/changes
    echo Build finished. The overview file is in %BUILDDIR%/changes.
    exit /b 0
)

if "%1" == "linkcheck" (
    %SPHINXBUILD% -b linkcheck . %BUILDDIR%/linkcheck
    echo Link check complete. Check the output in %BUILDDIR%/linkcheck/output.txt.
    exit /b 0
)

if "%1" == "doctest" (
    %SPHINXBUILD% -b doctest . %BUILDDIR%/doctest
    echo Doctests finished. Check the results in %BUILDDIR%/doctest/output.txt.
    exit /b 0
)

echo Unknown target: %1
exit /b 1