\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{spral}[2014/03/06 SPRAL Documentation Class]

% Load base
\LoadClass{article}

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

% headers and footerrs
\fancypagestyle{firststyle}
{
   \fancyhf{}
   %\lfoot{\texttt{module \MakeLowercase{\libraryname\_\packagename}}}
   \lfoot{\libraryurl}
   \cfoot{\thepage}
   \rfoot{\theversiondate}
   \renewcommand{\headrulewidth}{0pt}
}
\fancypagestyle{followstyle}
{
   \lhead{\bfseries \textcolor{stfcblue}{\libraryname} \texttt{\thepackagename}}
   %\lfoot{\texttt{module \MakeLowercase{\libraryname\_\packagename}}}
   \chead{Version \theversion}
   \lfoot{\libraryurl}
   \cfoot{\thepage}
   \rfoot{\theversiondate}
}


% Beginning of document hooks
\AtBeginDocument{
   % Set headers
   \pagestyle{followstyle}
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
\RequirePackage{tikz}
\newcommand{\funder}[2]{
   \hfill
   \begin{tikzpicture}
      \node (img) {\pgfimage[height=1cm]{#1}};
      \node[right,text width=3cm] at (img.east) {
         \footnotesize
         #2
      };
   \end{tikzpicture}
   \hfill
}

% New command definitions
\newcommand{\intentin}{\texttt{INTENT(IN)}}
\newcommand{\intentinout}{\texttt{INTENT(INOUT)}}
\newcommand{\intentout}{\texttt{INTENT(OUT)}}
\newcommand{\EqnRef}[1]{(\ref{#1})}

% Define our own \maketitle to avoid multiple call issues
\makeatletter
\newcommand{\thepkglang}{}
\newcommand{\pkglang}[1]{
   \renewcommand{\thepkglang}{#1}
}
\newcommand{\spralmaketitle}{
   \pgfimage[height=1cm]{stfc}
   \hfill
   {\Huge \bfseries \textcolor{stfcblue}{\libraryname}}
   \vspace{0.1cm}
   \textcolor{stfcgrey}{\hrule}
   \vspace{0.5mm}
   \hfill \textcolor{stfcgrey}{\Large \thepackagename\ v\theversion}
   \vspace{-0.3cm}

   \begin{center}
      {
         \LARGE \bfseries
         \texttt{\libraryname\_\thepackagename}\\
         \vspace{0.5em}
         \@title
      }

      \vspace*{0.5em}
      {\LARGE \thepkglang\ User Guide}
   \end{center}
   \begin{quote}
      \large
      \thepurpose
   \end{quote}

   \begin{flushright}
   %\noindent
   %{\bfseries Authors:}
   %\vspace{0.1cm}

   \noindent
   \@author
   \end{flushright}
   \vspace{-0.2cm}
   \textcolor{stfcgrey}{\hrule}
}
\makeatother

% Installation instruction need special maketitle alternative
\makeatletter
\newcommand{\makeinstalltitle}{
   \pgfimage[height=1cm]{stfc}
   \hfill
   {\Huge \bfseries \textcolor{stfcblue}{\libraryname}}
   \vspace{0.1cm}
   \textcolor{stfcgrey}{\hrule}
   \vspace{0.5cm}

   \begin{center}
      \LARGE \bfseries
      \@title
   \end{center}
   \begin{quote}
      \large
      \thepurpose
   \end{quote}

   \begin{flushright}
   \noindent
   \@author
   \end{flushright}
   \vspace{-0.2cm}
   \textcolor{stfcgrey}{\hrule}
}
\makeatother

\newenvironment{funders}{
   \vspace{1cm}
   \textcolor{stfcgrey}{\hrule}
   \vspace{0.1cm}
}{
   \vspace{0.1cm}
   \textcolor{stfcgrey}{\hrule}
}

% Avoid hyphenation
\hyphenpenalty=10000
%\sloppy

\endinput
