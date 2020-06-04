#!/usr/bin/env python3

"""
Create ARTS grammar by processing Jinja template and
install syntax highlighting extension for VSCode
"""

import os
from shutil import copytree, rmtree
from jinja2 import Environment, FileSystemLoader
import pyarts.workspace


def main():
    # Get names of ARTS workspace variables and methods
    arts_variables = pyarts.workspace.variables.workspace_variables.keys()
    arts_methods = pyarts.workspace.methods.workspace_methods.keys()

    # Add ARTS workspace variables and methods to language syntax configuration
    env = Environment(loader=FileSystemLoader("syntaxes"))
    template = env.get_template("arts.tmLanguage.in.json")
    with open("syntaxes/arts.tmLanguage.json", "w") as fp:
        fp.write(
            template.render(arts_methods=arts_methods, arts_variables=arts_variables)
        )

    # Install extension
    homedir = os.getenv("HOME")
    if homedir:
        extensiondir = os.path.join(
            homedir, ".vscode", "extensions", "arts-control-lang"
        )
        if extensiondir != os.path.abspath(os.curdir):
            print("Installing into " + extensiondir)
            rmtree(extensiondir)
            copytree(".", extensiondir)
    else:
        print("Unable to determine home directory location.")
        print(
            "Please copy arts-control-lang directory manually to "
            "~/.vscode/extensions/"
        )


if __name__ == "__main__":
    main()
