# GitHub basics
In order to contribute to this project you need to be
[registered on GitHub](github-join).

## Configure your local ``git`` client
GitHub matches individual commits and users using their email. To ensure a 
proper matching you should [configure your local ``git`` client][git-setup]
to use the email used to register on GitHub:
```bash
$ git config --global user.name "John Doe"
$ git config --global user.email johndoe@example.com
```

We highly recommend to [setup a SSH key][github-ssh] to simplify the 
authentication with the GitHub server
(``Profile > Settings > SSH and GPG keys > New SSH key``).

## Configure a new project

1. [Fork the project](#fork-the-project)
2. [Clone your fork](#clone-the-forked-project)
3. [Set up remotes](#set-up-remotes)
4. [Update your working copy](#update-your-working-copy)

### Fork the project
The first step is to [fork the project][github-fork] (click the ``Fork`` 
button on the project site). A fork is your own copy of the repository that is
stored server-side.

### Clone the forked project
Afterwards you can [clone your fork][atlassian-clone] to create a local 
working copy:
```bash
$ git clone https://github.com/YOUR_USERNAME/PROJECT.git
$ git clone git@github.com:YOUR_USERNAME/PROJECT.git  # With SSH key
```

### Set up remotes
The local ``git`` client is able to communicate with different remote 
repositories (remotes). It is convention to name the remote pointing to your
 fork ``origin`` (``git clone`` does that by default) and the remote 
 pointing to the original project repository ``upstream``.
 
You can [add the project URL][github-remote] to your local repository to be 
able to pull changes from it:
```bash
$ git remote add upstream https://github.com/atmtools/PROJECT.git
```

List all remotes to check the configuration of your ``git`` client:
```bash
$ git remote -v
origin https://github.com/YOUR_USERNAME/PROJECT.git (fetch)
origin https://github.com/YOUR_USERNAME/PROJECT.git (push)
upstream https://github.com/atmtools/PROJECT.git (fetch)
upstream https://github.com/atmtools/PROJECT.git (push)
```

### Update your working copy
Make sure to [pull changes][atlassian-pull] from the ``upstream`` ``master`` 
branch at regular intervals to keep track of changes done to the project.
We recommend to use the ``--rebase`` flag. This will
[rewind your commits and reapply them][atlassian-rebase] on top of the 
project's history to maintain a linear history.
```bash
$ git pull --rebase upstream master
```

[atlassian-rebase]: https://www.atlassian.com/git/tutorials/rewriting-history/git-rebase
[atlassian-clone]: https://www.atlassian.com/git/tutorials/setting-up-a-repository/git-clone
[atlassian-pull]: https://www.atlassian.com/git/tutorials/syncing/git-pull
[github-fork]: https://help.github.com/en/articles/fork-a-repo
[github-join]: https://github.com/join
[github-remote]: https://help.github.com/en/articles/adding-a-remote
[github-ssh]: https://help.github.com/en/articles/connecting-to-github-with-ssh
[git-setup]: https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup
