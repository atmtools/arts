# -*- coding: utf-8 -*-

"""This module contains custom roles to use in Sphinx.
"""
from docutils import nodes


def setup(app):
    """Install the extension.

    Parameters:
        app: Sphinx application context.
    """
    app.add_role('arts', arts_docserver_role)


def arts_docserver_role(name, rawtext, text, lineno, inliner, options=None,
                        content=None):
    """Create a link to ARTS docserver.

    Parameters:
        name (str): The role name used in the document.
        rawtext (str): The entire markup snippet, with role.
        text (str): The text marked with the role.
        lineno (str): The line number where rawtext appears in the input.
        inliner (str): The inliner instance that called us.
        options (dict): Directive options for customization.
        content (list): The directive content for customization.

    Returns:
     list, list: Nodes to insert into the document, System messages.
    """
    if content is None:
        content = []

    if options is None:
        options = {}

    url = 'https://atmtools.github.io/arts-docs-master/docserver/all/{}'.format(text)
    node = nodes.reference(rawtext, text, refuri=url, **options)

    return [node], []
