#!/usr/bin/env python

# Create ARTS grammar by processing Jinja template

from jinja2 import Environment, FileSystemLoader
import pyarts.workspace

def main():
    # Get names of ARTS workspace variables and methods
    arts_variables = pyarts.workspace.variables.workspace_variables.keys()
    arts_methods = pyarts.workspace.methods.workspace_methods.keys()
    
    # Add ARTS workspace variables and methods to language syntax configuration
    env = Environment(loader=FileSystemLoader('syntaxes'))
    template = env.get_template("arts.tmLanguage.in.json")
    with open("syntaxes/arts.tmLanguage.json","w") as fp:
        fp.write(template.render(arts_methods=arts_methods,arts_variables=arts_variables))

if __name__ == "__main__":
    main()
