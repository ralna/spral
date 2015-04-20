\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{spralweb}[2014/03/06 SPRAL Documentation Class]

% Load base
\LoadClass{book}

% Library name definition
\RequirePackage{url}
\newcommand{\libraryname}{SPRAL}
\newcommand{\libraryurl}{\url{http://www.numerical.rl.ac.uk/spral}}
\newcommand{\spral}{\texttt{\libraryname}}

% Colors
\RequirePackage{pgf}
\definecolor{stfcblue}{RGB}{21,42,75}
\definecolor{stfcgrey}{RGB}{152,155,151}

% Required packages
\RequirePackage{fancyhdr}
\RequirePackage[top=2cm,left=2.5cm,right=2.5cm,bottom=2cm,headheight=12pt]{geometry}
\RequirePackage{verbatim}
\RequirePackage{algorithmic}

% Metadata commands
\newcommand{\thepackagename}{}
\newcommand{\packagename}[1]{
   \renewcommand{\thepackagename}{#1}
}
\newcommand{\theversion}{}
\newcommand{\version}[1]{
   \renewcommand{\theversion}{#1}
}
\newcommand{\theversiondate}{}
\newcommand{\versiondate}[1]{
   \renewcommand{\theversiondate}{#1}
}
\newcommand{\thepurpose}{}
\newcommand{\purpose}[1]{
   \renewcommand{\thepurpose}{#1}
}
\newcommand{\thepkglang}{}
\newcommand{\pkglang}[1]{
   \renewcommand{\thepkglang}{#1}
}

% Package information headers
\newcommand{\authors}{
}
\newcommand{\versionhistory}{
   \vspace{0.6cm}
   \noindent
   {\large \bfseries Major Version History}
}

% Support for specifying funders
\newenvironment{funders}{
   \HCode{\Hnewline<div class="funders">\Hnewline}
}{
   \HCode{\Hnewline</div> <!-- funders -->\Hnewline}
}
\newcommand{\funder}[2]{
   \HCode{<div class="funder">\Hnewline}
   \HCode{<img src="#1.png">\Hnewline<p>#2</p>\Hnewline}
   \HCode{</div> <!-- funder -->\Hnewline}
}

% New command definitions
\newcommand{\intentin}{\texttt{INTENT(IN)}}
\newcommand{\intentinout}{\texttt{INTENT(INOUT)}}
\newcommand{\intentout}{\texttt{INTENT(OUT)}}
\newcommand{\EqnRef}[1]{(\ref{#1})}

% Define our own \maketitle to avoid multiple call issues
\makeatletter
\newcommand{\spralmaketitle}{
   \HCode{<div class="maketitle">\Hnewline}
   \HCode{\Hnewline<h2 class="pkgName">}\libraryname\_\thepackagename\HCode{<span class="pkgVersion">}v\theversion\HCode{</span></h2>\Hnewline}
   \HCode{<h2 class="titleHead">}
      \@title
   \HCode{</h2>\Hnewline}
   \HCode{<span class="pkgLang">}
      \thepkglang\ User Guide
   \HCode{</span> <br />\Hnewline}
   \HCode{<div class="pkgPurpose">\Hnewline}
      \thepurpose
   \HCode{\Hnewline</div> <!-- pkgPurpose -->\Hnewline}
   \HCode{<div class="author">}
      \@author
   \HCode{</div>\Hnewline}
   \HCode{</div>\Hnewline}
}
\makeatother

% Installation instruction need special maketitle alternative
\makeatletter
\newcommand{\makeinstalltitle}{
   \HCode{<div class="maketitle">\Hnewline}
   \HCode{<h2 class="titleHead">}
      \@title
   \HCode{</h2>\Hnewline}
   \HCode{<div class="pkgPurpose">\Hnewline}
      \thepurpose
   \HCode{\Hnewline</div> <!-- pkgPurpose -->\Hnewline}
   \HCode{<div class="author">}
      \@author
   \HCode{</div>\Hnewline}
   \HCode{</div>\Hnewline}
}
\makeatother

% Avoid hyphenation
\hyphenpenalty=10000
%\sloppy

\endinput
