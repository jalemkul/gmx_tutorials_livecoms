# Preparing LaTex articles for the Living Journal of Computational Molecular Sciences ([LiveCoMS](http://www.livecomsjournal.org/)) using GitHub

## LaTeX and GitHub

LaTeX is a powerful typesetting program that is very commonly used in the mathematics and physics community, and increasingly in the computational sciences.
We have provided a LaTeX template for use with LiveCoMS articles in order to make it easy to have a pleasing, consistent visual style.

We request that articles developed for LiveCoMS do so by creating a GitHub repository for the paper.
GitHub was created to be able to easily manage code versioning, merging, and branching for software development.
We feel that for collaborative documents like the ones in LiveCoMS, which encourage feedback from the community, and not just for authors, this software model is a very good fit.
See a more on-depth discussion of this in [brief explanation of "paper writing as code development"](http://https://github.com/livecomsjournal/journal_information/blob/master/editorial_material/PAPER_CODE.md) for this model we are using for LiveCoMS.

### Setting up and Managing a GitHub repository

GitHub has a substantial set of online tutorials and answered questions.  [It is very easy to set up a repository in GitHub](https://help.github.com/articles/create-a-repo/).

Some particular recommendations we have when you set up your GitHub repository.
- Use the issue tracker.  The [issue tracker model](https://guides.github.com/features/issues/) is a powerful way to manage contributions from a number of different authors and contributors, and allow community engagement with your work.
- [Protect the master branch](https://help.github.com/articles/about-protected-branches/); if people are changing the master branch directly, it will be much harder to track revisions and make sure changes are edited before being incorporated.
- We recommend adding changes to the repository [through pull requests](https://help.github.com/articles/proposing-changes-to-your-work-with-pull-requests/). In this way, other people can easily see what changes are being requested and made. Pull requests can be made from a branch of your repository (if from one of your co-authors) or from a fork of your repository (if from an outside individual who would like to contribute) and are handled similarly in either case -- the repository owners review and decide whether or not to acccept the proposed changes.

### Versioning and dates

All versions of your article should be dated; our templates automatically add a date to each compilation of the PDF, and when you make changes to your content, you should include an updated version of the PDF in your repository.

**Versioning of articles occurs via two mechanisms which should work together**:
- LiveCoMS templates use titles which have a version number explicitly stated in the title. Change this version whenever you make a significant change, using a numbering scheme you agree on with your co-authors and, ideally, document here.
- GitHub has a mechanism to archive specific versions of code (in this case, your paper!) as specific "releases". Whenever you increment the version number in your title, you should also [create a new release on GitHub](https://help.github.com/articles/creating-releases/) to archive it and concisely describe any associated changes.

It is also possible to easily set up tools like Zenodo to automatically generate a new DOI for each release of a project, if you so desire.

#### Versioning and LiveCoMS

Please note that versioning on GitHub is your choice and is separate from the process of publishing new peer-reviewed versions in LiveCoMS.
Each version you publish in LiveCoMS must have a unique version number, but these version numbers are up to you, and can track with what is used on your GitHub repository.
For example, you might choose to give your first published version on LiveCoMS the version number 1.0, but if in two years it makes sense to do a new peer-reviewed version, this could hypothetically be version 2.4, if that's what your GitHub repository is at.

#### Change logs

We recommend you keep track of even minor changes in a running log near the bottom of the base `README.md` file in your repository, and any time enough minor changes accumulate, create a new version number and release to clearly indicate that these have been incorporated and update your revision log.
This will allow readers to easily track progress of your document over time without having to navigate the full commit history on GitHub.

#### Example scenarios for authors

Perhaps you've just come out with a new Best Practices article for LiveCoMS, and several community members are very interested and contact you via your issue tracker with ideas for other material which could be included, and another contacts you by e-mail.
We strongly recommend funneling as much discusion as possible to your GitHub issue tracker since it will allow interested individuals to discuss publicly and openly, which e-mail does not.
So, first, we suggest referring the e-mail contact to GitHub (if he or she is interested).
For individuals who suggest including new material, if you are receptive, you have several options:
- Write material like what they suggest and incorporate it (acknowledging them as appropriate)
- See if they are interested in actually providing draft material for you to incorporate; this could come in as a pull request (where they request to change your document and you review and propose modifications and accept it if you're happy with it) or they could just provide you with content. The former is preferable since GitHub will then formally track their contribution. One issue to consider in this case is whether their contributions will merit addition as an author. You should work this out with them before they provide significant draft content. Note that you CAN add authors to your document, though authorship will not change at LiveCoMS until a new peer reviewed version is published.

In either case you would close the relevant GitHub issue once addressed.

#### Example scenarios for readers

Perhaps you've just found a LiveCoMS article which interests you greatly, and you see some important references it missed or have an idea for other material which really ought to be included.
We recommend that your first step be to raise an issue on the repository's issue tracker, discussing your observations and asking the authors' opinion about how to handle it.
You should be sure to indicate the nature of your interest.
For example, you might say:
- I just wanted to make sure you were aware of this in case you would like to discuss it
- I believe you should write new content which discusses and/or cites this for the following reasons
- I'm very interested in this work and would be happy to contribute to it by adding material along these lines if you are interested, or if not, I could give you some suggestions on relevant references if you would like to investigate

Basically, try to give authors some sense of how interested you are, how strongly you feel about the issue, and whether you are willing to get involved or just providing feedback.
Be aware that some authors may be delighted to add new authors and contributions to a paper, whereas others may want to keep the authors list short and tightly curate everything.
DO make sure that, if you head towards actually writing any content for inclusion, everyone's expectations about authorship are clear in advance.

If you do end up writing something, you will likely need to make a fork of the repository, make your changes on that fork, and then submit a pull request from your fork in order to change the repository.
If you head this route, be sure you have read the [LiveCoMS authors guide](https://livecomsjournal.github.io/authors/), as you want to make sure your changes continue to conform to what LiveCoMS requires.

### Some LaTeX formatting tips

- Keep your paper in one sentence per line.  This makes it easy to see which sentences have changed when you look at the diff between one version and the next.  The GitHub diff tools operate on a per-line basis. If there are no line breaks between sentences, then other authors will see only that something somewhere in the paragraph has changed.
