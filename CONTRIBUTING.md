# How to contribute
Thank you for considering to contribute to our project! :tada: :+1:

In order to contribute you need to be
[registered on GitHub](https://github.com/join) and familiar with the
concepts of [forking and communicating](GITHUB-INTRO.md) with GitHub.

## Checklist
This checklist is closely based on a
[blog post](http://blog.davidecoppola.com/2016/11/howto-contribute-to-open-source-project-on-github/)
by Davide Coppola. It also serves as a table of contents.

1. [Create a branch](#create-a-branch)
2. [Work on your contribution](#work-on-your-contribution)
3. [Pull request](#pull-request)
4. [Code review](#code-review)
5. [Clean up](#clean-up)

## Create a branch
Before your first contribution, make sure to
[fork the project](GITHUB-INTRO.md#configure-a-new-project).

Create a local branch to collect changes related to a new feature or bug fix. 
Branches help to organize different developments within the project and 
facilitate the code review process.

You can do that with the following git command:
```bash
$ git checkout -b BRANCH_NAME
```

This will [create][atlassian-branch] and [checkout][atlassian-checkout]
a new branch ``BRANCH_NAME`` in your local repository.
It is best practice to use a name that clearly describes the purpose of the 
branch.

You can check which branch is active (``*``) using the ``git`` client:
```bash
$ git branch
  master
* BRANCH_NAME
```

## Work on your contribution
Now you can start with the development of your new feature (or bug fix).
During the work you can use ``git``'s commit and push mechanism to save and
track your changes to your personal fork:
```bash
$ git push origin BRANCH_NAME
```

## Pull request
After pushing your changes to your fork navigate to the GitHub page of your
fork and [create a Pull request][github-pr]. Add the needed information to the
web form and submit your request.

## Code review
The developer team will review your changes and decide whether to accept your
changes. This process might include some discussion or even further changes to
the code (this is the reason why [branches](#create-a-branch) are important).

## Clean up
After your contribution has been merged to the main project (or rejected) you
can delete your development branch:
```bash
$ git branch -D BRANCH_NAME  # local
$ git push origin --delete BRANCH_NAME  # server-side
```

[atlassian-branch]: https://www.atlassian.com/git/tutorials/using-branches/
[atlassian-checkout]: https://www.atlassian.com/git/tutorials/using-branches/git-checkout
[github-pr]: https://help.github.com/en/articles/creating-a-pull-request
