;;; gmic-mode.el --- Major mode for editing G'MIC scripts -*- lexical-binding: t -*-

;; Author: Generated for G'MIC scripting
;; Version: 0.1
;; Keywords: languages, G'MIC, image processing
;; URL: https://gmic.eu

;;; Commentary:
;; A major mode for editing G'MIC script files (.gmic).
;; Features:
;;  - Syntax highlighting (font-lock) adapted to G'MIC syntax
;;  - Indentation based on G'MIC block structures
;;  - Comment handling with '#'
;;  - Navigation between command definitions (imenu)
;;  - Optional: run current script with M-x gmic-run

;;; Code:


;;;; -------------------------------------------------------------------------
;;;; Customization

(defgroup gmic nil
  "Major mode for editing G'MIC scripts."
  :group 'languages
  :prefix "gmic-")

(defcustom gmic-indent-offset 2
  "Number of spaces per indentation level in G'MIC scripts."
  :type 'integer
  :group 'gmic)

(defcustom gmic-executable "gmic"
  "Path to the G'MIC executable used by `gmic-run'."
  :type 'string
  :group 'gmic)

;;;; -------------------------------------------------------------------------
;;;; Syntax table

(defvar gmic-mode-syntax-table
  (let ((st (make-syntax-table)))
    ;; '#' starts a line comment
    (modify-syntax-entry ?# "<" st)
    (modify-syntax-entry ?\n ">" st)
    ;; String delimiters
    (modify-syntax-entry ?\" "\"" st)
    ;; Braces and brackets as paired delimiters
    (modify-syntax-entry ?{ "(}" st)
    (modify-syntax-entry ?} "){" st)
    (modify-syntax-entry ?\[ "(]" st)
    (modify-syntax-entry ?\] ")[" st)
    ;; '$' is part of variable names
    (modify-syntax-entry ?$ "_" st)
    ;; Some operators
    (modify-syntax-entry ?% "." st)
    (modify-syntax-entry ?& "." st)
    (modify-syntax-entry ?| "." st)
    st)
  "Syntax table for `gmic-mode'.")

;;;; -------------------------------------------------------------------------
;;;; Faces

(defface gmic-command-name-face
  '((t :inherit font-lock-function-name-face :weight bold))
  "Face used for G'MIC command definition names."
  :group 'gmic)

;;;; -------------------------------------------------------------------------
;;;; Font-lock keywords

(defvar gmic--block-open-keywords
  '("repeat" "for" "foreach" "do" "if" "elif" "else" "local" "l")
  "G'MIC keywords that open an indented block.")

(defvar gmic--block-close-keywords
  '("done" "while" "fi" "elif" "else")
  "G'MIC keywords that close an indented block.")

(defvar gmic--control-keywords
  '("if" "elif" "else" "fi"
    "repeat" "for" "foreach" "done"
    "do" "while"
    "local" "l"
    "skip" "break" "continue"
    "return" "error" "warning"
    "pass")
  "G'MIC control flow keywords.")

(defvar gmic--builtin-commands
  ;; Native G'MIC commands, hard-coded in C++ — stable list.
  ;; Purely symbolic operators (!=, %, &, *, +, -, /, etc.)
  ;; are handled separately via their own font-lock rule.
  '("abs" "abscut" "acos" "acosh" "add" "add3d" "and"
    "append" "asin" "asinh" "atan" "atan2" "atanh"
    "b" "bilateral" "blur" "boxfilter" "break" "bsl" "bsr"
    "c" "camera" "check" "check3d" "command" "continue"
    "convolve" "correlate" "cos" "cosh" "createdir" "crop"
    "cumulate" "cursor" "cut"
    "debug" "delete" "denoise" "deriche" "dilate" "discard"
    "displacement" "distance" "div" "div3d" "do" "done"
    "e" "echo" "eigen" "elif" "ellipse" "else" "endian"
    "eq" "equalize" "erf" "erode" "error" "eval" "exec" "exp"
    "f" "fft" "fi" "files" "fill" "flood" "for" "foreach"
    "g" "ge" "gt" "guided"
    "histogram"
    "i" "if" "ifft" "image" "index" "inpaint" "input"
    "isoline3d" "isosurface3d"
    "j" "j3d"
    "k" "keep"
    "l" "l3d" "label" "le" "light3d" "line" "local" "log" "log2" "log10" "lt"
    "m" "map" "matchpatch" "max" "maxabs" "mdiv" "median"
    "min" "minabs" "mirror" "mmul" "mod" "move" "mproj"
    "mul" "mul3d" "mv"
    "n" "name" "named" "neq" "network" "nm" "nmd" "noarg"
    "noise" "normalize"
    "o" "object3d" "onfail" "or" "output"
    "parallel" "pass" "permute" "point" "polygon" "pow" "progress"
    "q" "qr" "quit"
    "r" "r3d" "rand" "remove" "repeat" "resize" "return"
    "reverse" "rm" "rol" "ror" "rotate" "rotate3d" "round" "rv"
    "s" "screen" "serialize" "set" "sh" "shared" "shift"
    "sign" "sin" "sinc" "sinh" "skip" "smooth" "solve" "sort"
    "split" "sqr" "sqrt" "srand" "status" "store"
    "streamline3d" "sub" "sub3d" "svd"
    "t" "tan" "tanh" "text"
    "u" "uncommand" "unroll" "unserialize"
    "v" "vanvliet" "verbose"
    "w" "wait" "warn" "warp" "watershed" "while" "window"
    "x" "xor"
    "y"
    "z")
  "G'MIC native built-in commands, hard-coded in C++.")


(defconst gmic-font-lock-keywords
  (list
   ;; 1. Command definitions: "name :" or "name(args) :"
   ;;    ':' must not be followed by '=' (otherwise it is a var:= assignment).
   ;;    Use :[^=] or ':' at end of line (no lookahead in Emacs regexp).
   '("^\\s-*\\([a-zA-Z_][a-zA-Z0-9_]*\\)\\s-*\\(([^)]*)\\s-*\\)?:\\([^=]\\|$\\)"
     (1 'gmic-command-name-face))

   ;; 2. Control flow keywords (standalone on a line or with arguments)
   (cons (concat "\\_<"
                 (regexp-opt gmic--control-keywords)
                 "\\_>")
         'font-lock-keyword-face)

   ;; 3. Variable references: $var, ${var}, $>, $<, $!, $#, $$
   '("\\(\\$\\(?:{[^}]*}\\|[a-zA-Z_][a-zA-Z0-9_]*\\|[><#!$0-9]\\)\\)"
     (1 font-lock-variable-name-face))

   ;; 4. Image selectors: [n], [-n], [name], [#$var], [^n]
   ;;    (inside brackets after a command)
   '("\\[\\(\\^?-?[0-9]+\\|\\^?[a-zA-Z_][a-zA-Z0-9_]*\\|#[$a-zA-Z0-9_]*\\|[,0-9 -]*\\)\\]"
     (0 font-lock-constant-face))

   ;; 5. Assignment operators := and =>
   '("[:=]=" . font-lock-builtin-face)
   '("=>" . font-lock-builtin-face)

   ;; 6. Command prefix operators: +cmd, -cmd, *cmd at word boundary
   '("\\([+\\-\\*/]\\)\\([a-zA-Z_][a-zA-Z0-9_]*\\)"
     (1 font-lock-preprocessor-face)
     (2 font-lock-builtin-face))

   ;; 7. Native built-in commands — anywhere on the line, at word boundary
   (cons (concat "\\_<" (regexp-opt gmic--builtin-commands t) "\\_>")
         'font-lock-builtin-face)

   ;; 8. Symbolic native operators: +3d -3d *3d /3d m* m/ => != == <= >= << >>
   '("\\(?:\\+3d\\|-3d\\|\\*3d\\|/3d\\|m\\*\\|m/\\|!=\\|==\\|<=\\|>=\\|<<\\|>>\\|=>\\)"
     . font-lock-builtin-face)

   ;; 8. Numeric literals (integers and floats, including degree notation)
   '("\\b\\([0-9]+\\.?[0-9]*\\(?:e[+-]?[0-9]+\\)?°?\\)\\b"
     (1 font-lock-constant-face))

   ;; 9. Shebang line
   '("^#!.*$" . font-lock-comment-face)

   ;; 10. Special constants
   '("\\_<\\(pi\\|inf\\|nan\\|true\\|false\\)\\_>" . font-lock-constant-face))
  "Font-lock keywords for `gmic-mode'.")

;;;; -------------------------------------------------------------------------
;;;; Indentation — counting opening/closing keywords per line

(defvar gmic--indent-open-re
  ;; Only keywords count as openers.
  ;; '{' is never counted: it is a syntactic shortcut for the keyword
  ;; that precedes it, which is already counted.
  (concat "\\_<"
          (regexp-opt '("repeat" "for" "foreach" "do"
                        "if" "elif" "else"
                        "local" "l")
                      t)
          "\\_>")
  "Regexp matching block-opening keywords in G'MIC.")

(defvar gmic--indent-close-re
  ;; Isolated '}' (preceded by whitespace or start of line, followed by whitespace or end)
  ;; plus closing keywords.
  (concat "\\(?:\\(?:^\\|\\s-\\)}\\(?:\\s-\\|$\\)\\|\\_<"
          (regexp-opt '("done" "fi" "elif" "else" "while") t)
          "\\_>\\)")
  "Regexp matching block-closing keywords or isolated '}' in G'MIC.")

(defun gmic--comment-start-p (str i)
  "Return t if '#' at position I in STR starts a comment.
A '#' is a comment only when at line start (after whitespace) or
preceded by a space or tab.  A '#' immediately following a non-blank
character is a G'MIC image selector (e.g. '#$var', '#0')."
  (or (= i 0)
      (let ((prev (aref str (1- i))))
        (or (= prev ?\s) (= prev ?\t)))))

(defun gmic--strip-comments-and-strings (line-str)
  "Return LINE-STR with comments and string contents removed.
Keeps the quote delimiters themselves to preserve quote counting,
but removes everything between them.
A '#' is only a comment when preceded by a space/tab or at line start."
  (let ((result "")
        (in-string nil)
        (i 0)
        (len (length line-str)))
    (while (< i len)
      (let ((ch (aref line-str i)))
        (cond
         ((= ch ?\")
          (setq in-string (not in-string))
          (setq result (concat result "\"")))
         ((and (not in-string)
               (= ch ?#)
               (gmic--comment-start-p line-str i))
          (setq i len))               ; stop — rest is comment
         ((not in-string)
          (setq result (concat result (string ch))))))
      (setq i (1+ i)))
    result))

(defun gmic--strip-comments (line-str)
  "Return LINE-STR with any trailing G'MIC comment removed.
String contents are preserved (unlike `gmic--strip-comments-and-strings').
A '#' is only a comment when preceded by a space/tab or at line start."
  (let ((result "")
        (in-string nil)
        (i 0)
        (len (length line-str)))
    (while (< i len)
      (let ((ch (aref line-str i)))
        (cond
         ((= ch ?\")
          (setq in-string (not in-string))
          (setq result (concat result "\"")))
         ((and (not in-string)
               (= ch ?#)
               (gmic--comment-start-p line-str i))
          (setq i len))               ; stop — rest is comment
         (t
          (setq result (concat result (string ch))))))
      (setq i (1+ i)))
    result))


(defun gmic--count-occurrences (re str)
  "Count non-overlapping matches of RE in STR."
  (let ((count 0)
        (start 0))
    (while (string-match re str start)
      (setq count (1+ count)
            start (match-end 0)))
    count))

(defvar gmic--command-def-re
  "^\\s-*[a-zA-Z_][a-zA-Z0-9_]*\\s-*\\(?:([^)]*)\\s-*\\)?:\\([^=]\\|$\\)"
  "Regexp matching a G'MIC command definition line (e.g. 'my_cmd :').
The ':' must not be followed by '=' to avoid matching 'var:=value'.")

(defun gmic--line-is-command-def-p (line-str)
  "Return t if LINE-STR is a G'MIC command definition."
  (string-match-p gmic--command-def-re line-str))

(defun gmic--line-is-cli-gui-comment-p (line-str)
  "Return t if LINE-STR is a G'MIC special documentation comment.
Lines starting with `#@cli' or `#@gui' (possibly preceded by whitespace)
have a special status in G'MIC and must always be indented at column 0."
  (string-match-p "^[[:space:]]*#@\\(?:cli\\|gui\\)" line-str))

(defun gmic--line-delta (line-str)
  "Return the net indentation delta produced by LINE-STR.
Positive means the *next* line should be indented further.
Closing keywords (fi, done, while, elif, else, isolated '}') are NOT
counted here: their effect on the current line's position is handled by
`gmic--line-leading-close-delta', and they do not affect the next line.
Only opening keywords produce a positive delta.
A command definition line counts as net +1 (it opens a command body)."
  (let* ((stripped (gmic--strip-comments-and-strings line-str))
         (open-re (concat "\\_<"
                          (regexp-opt '("repeat" "for" "foreach" "do"
                                        "if" "elif" "else"
                                        "local" "l") t)
                          "\\_>"))
         (opens  (gmic--count-occurrences open-re stripped))
         ;; Closers on the same line as an opener cancel that opener.
         ;; e.g. "repeat 6 { ... done" on one line = delta 0.
         ;; Only count closers that appear AFTER the first opener on the line.
         (close-re (concat "\\(?:\\(?:^\\|\\s-\\)}\\(?:\\s-\\|$\\)\\|\\_<"
                           (regexp-opt '("done" "fi" "while") t)
                           "\\_>\\)"))
         ;; Closers appearing AFTER the first opener on the line
         (first-open-pos (if (string-match open-re stripped)
                             (match-beginning 0)
                           nil))
         (trailing-closes
          (if first-open-pos
              (gmic--count-occurrences close-re
                                       (substring stripped first-open-pos))
            0))
         (cmd-def (if (gmic--line-is-command-def-p line-str) 1 0)))
    (* (+ (- opens trailing-closes) cmd-def) gmic-indent-offset)))

(defun gmic--line-leading-close-delta (line-str)
  "Return the dedent to apply to LINE-STR itself.
When a line starts with closing keywords or isolated '}', pull it back
by one level per leading closer.
'elif' and 'else' dedent by exactly one level (back to the 'if' level)
but are not counted via the general closer mechanism."
  (let* ((stripped (gmic--strip-comments-and-strings line-str))
         ;; elif/else: always -1 level relative to the previous line
         (leading-elif-else
          (if (string-match-p
               (concat "^\\s-*\\_<\\(?:elif\\|else\\)\\_>")
               stripped)
              1 0))
         ;; Other closers: done, fi, while, isolated '}' — before the first opener
         (first-open
          (let ((pos (and (string-match
                           (concat "\\_<"
                                   (regexp-opt '("repeat" "for" "foreach" "do"
                                                 "if" "local" "l") t)
                                   "\\_>")
                           stripped)
                          (match-beginning 0))))
            (or pos (length stripped))))
         (substr (substring stripped 0 first-open))
         (other-close-re
          (concat "\\(?:\\(?:^\\|\\s-\\)}\\(?:\\s-\\|$\\)\\|\\_<"
                  (regexp-opt '("done" "fi" "while") t)
                  "\\_>\\)"))
         (leading-others (gmic--count-occurrences other-close-re substr)))
    (* (+ leading-elif-else leading-others) gmic-indent-offset)))

(defun gmic--count-unescaped-quotes (str)
  "Count double-quote characters in STR, stopping at an unquoted '#' comment."
  (let ((count 0)
        (in-string nil)
        (i 0)
        (len (length str)))
    (while (< i len)
      (let ((ch (aref str i)))
        (cond
         ((= ch ?\")
          (setq count (1+ count)
                in-string (not in-string)))
         ((and (not in-string) (= ch ?#))
          (setq i len)))) ; end of line: comment
      (setq i (1+ i)))
    count))

(defun gmic--in-multiline-string-p ()
  "Return t if point is currently inside a multiline G'MIC string.
Counts all double-quote characters from the start of the current
command definition (or buffer start) up to the previous line.
An odd total means we are inside an open string."
  (let ((total-quotes 0))
    (save-excursion
      (beginning-of-line)
      ;; Walk back to the enclosing command definition or buffer start
      (let ((limit (save-excursion
                     (if (re-search-backward gmic--command-def-re nil t)
                         (line-beginning-position)
                       (point-min)))))
        (while (> (point) limit)
          (forward-line -1)
          (let ((line (gmic--current-line-string)))
            (unless (string-match-p "^\\s-*\\(?:#.*\\)?$" line)
              (setq total-quotes
                    (+ total-quotes
                       (gmic--count-unescaped-quotes line))))))))
    (= (% total-quotes 2) 1)))

(defun gmic--current-line-string ()
  "Return the content of the current line as a string."
  (buffer-substring-no-properties
   (line-beginning-position) (line-end-position)))

(defun gmic--line-indentation (line-str)
  "Return the number of leading spaces in LINE-STR."
  (let ((i 0))
    (while (and (< i (length line-str))
                (= (aref line-str i) ?\s))
      (setq i (1+ i)))
    i))

(defun gmic--paren-delta (line-str)
  "Return the net paren/bracket delta for LINE-STR.
Counts '(' and '[' as openers, ')' and ']' as closers.
Characters inside nested strings (double-quoted) are ignored.
A '#' is only treated as a comment when preceded by a space/tab or at
line start; '#$var' and '#N' are G'MIC image selectors, not comments."
  (let ((delta 0)
        (in-string nil)
        (i 0)
        (len (length line-str)))
    (while (< i len)
      (let ((ch (aref line-str i)))
        (cond
         ((= ch ?\")
          (setq in-string (not in-string)))
         ((not in-string)
          (cond
           ((or (= ch ?\() (= ch ?\[)) (setq delta (1+ delta)))
           ((or (= ch ?\)) (= ch ?\])) (setq delta (1- delta)))
           ((= ch ?#)
            (when (gmic--comment-start-p line-str i)
              (setq i len)))))))  ; stop — rest is comment
      (setq i (1+ i)))
    (* delta gmic-indent-offset)))

(defun gmic--paren-leading-close-delta (line-str)
  "Return the dedent for LINE-STR itself due to leading ')' or ']'.
Counts consecutive closing parens/brackets (possibly separated by
spaces) that appear before any opener on the line."
  (let ((delta 0)
        (i 0)
        (len (length line-str)))
    ;; Skip leading whitespace
    (while (and (< i len)
                (let ((ch (aref line-str i)))
                  (or (= ch ?\s) (= ch ?\t))))
      (setq i (1+ i)))
    ;; Count leading closers — stop at anything that is not ) ] or space
    (while (< i len)
      (let ((ch (aref line-str i)))
        (cond
         ((or (= ch ?\s) (= ch ?\t))
          (setq i (1+ i)))
         ((or (= ch ?\)) (= ch ?\]))
          (setq delta (1+ delta))
          (setq i (1+ i)))
         (t
          (setq i len)))))
    (* delta gmic-indent-offset)))

(defun gmic--find-matching-open-indent (line-str)
  "For a LINE-STR starting with ')' or ']', find the indentation of the
line that opened the corresponding paren/bracket.
Returns nil if not inside a multiline string or line doesn't start with a closer."
  (let* ((stripped (string-trim-left line-str))
         (first-char (and (> (length stripped) 0) (aref stripped 0))))
    (when (and first-char (or (= first-char ?\)) (= first-char ?\])))
      ;; Count how many leading closers we need to match
      (let ((need 0)
            (i 0)
            (len (length stripped)))
        (while (and (< i len)
                    (let ((ch (aref stripped i)))
                      (or (= ch ?\)) (= ch ?\]) (= ch ?\s) (= ch ?\t))))
          (let ((ch (aref stripped i)))
            (when (or (= ch ?\)) (= ch ?\]))
              (setq need (1+ need))))
          (setq i (1+ i)))
        (save-excursion
          (beginning-of-line)
          (let ((balance need)
                (result nil))
            (while (and (> balance 0) (not (bobp)))
              (forward-line -1)
              (let ((prev (gmic--current-line-string)))
                (unless (string-match-p "^\\s-*\\(?:#.*\\)?$" prev)
                  (let* ((d (/ (gmic--paren-delta prev) gmic-indent-offset))
                         (new-balance (- balance d)))
                    (when (<= new-balance 0)
                      ;; The closing paren aligns with the indentation of
                      ;; the line that opened the corresponding block.
                      (setq result (gmic--line-indentation prev)))
                    (setq balance new-balance)))))
            result))))))



(defun gmic--line-opens-string-p (line-str)
  "Return t if LINE-STR opens a multiline string (odd number of quotes)."
  (= (% (gmic--count-unescaped-quotes line-str) 2) 1))

(defun gmic--line-is-closing-quote-p (line-str)
  "Return t if LINE-STR consists solely of a closing double-quote
\(possibly surrounded by whitespace).  Such a line ends a multiline
string and should be indented at the same level as the line that
opened it, not at the indentation level of the string contents."
  (string-match-p "^[[:space:]]*\"[[:space:]]*$" line-str))

(defun gmic--find-string-opener-indent ()
  "Scan backward from the current line to find the line that opened the
current multiline string, and return its indentation column.
Returns 0 if the opener cannot be found."
  (save-excursion
    (beginning-of-line)
    (let ((quotes 0)
          (result 0)
          (found nil))
      ;; Walk backward; accumulate quote counts.  When the running total
      ;; becomes odd we have found the line that contains the unmatched
      ;; opening quote.
      (while (and (not found) (not (bobp)))
        (forward-line -1)
        (let ((line (gmic--current-line-string)))
          (unless (string-match-p "^\\s-*\\(?:#.*\\)?$" line)
            (setq quotes (+ quotes (gmic--count-unescaped-quotes line)))
            (when (= (% quotes 2) 1)
              (setq result (gmic--line-indentation line))
              (setq found t)))))
      result)))


(defun gmic-indent-line ()
  "Indent the current line for `gmic-mode'.

Outside multiline strings:
  1. Find the previous non-empty line that is not inside a multiline string.
  2. Start from that line's indentation + keyword delta.
  3. Subtract the leading-close contribution of the current line.

Inside multiline strings (math expressions):
  1. If the line starts with ')' or ']', locate the matching opener and
     align with its indentation.  This check is done before the backward
     scan so that point is still on the current line when the search runs.
  2. Otherwise find the previous non-empty line, start from its indentation
     + paren/bracket delta, and subtract leading closers.

In both cases: clamp to zero.
Special cases forced to column 0 regardless of context:
  - G'MIC command definitions (name :)
  - `#@cli' and `#@gui' documentation comment lines."
  (interactive)
  (let* ((in-string (gmic--in-multiline-string-p))
         (current-line (gmic--current-line-string))
         (indent 0))
    ;; Special case: a line containing only a closing `"` ends a multiline
    ;; string.  It must be indented at the same column as the line that
    ;; opened the string, not at the indentation of the string body.
    (if (and in-string (gmic--line-is-closing-quote-p current-line))
        (indent-line-to (gmic--find-string-opener-indent))
    (save-excursion
      (beginning-of-line)
      ;; Inside a multiline string a line beginning with ) or ] must align
      ;; with the line that opened the matching paren/bracket.  We resolve
      ;; this HERE, while point is still on the current line, because
      ;; gmic--find-matching-open-indent scans backward from (point) and
      ;; would produce wrong results if called after forward-line -1 below.
      (let ((matching (when in-string
                        (gmic--find-matching-open-indent current-line))))
        (if matching
            (setq indent matching)
          ;; General case: scan backward for the previous non-empty line.
          (let ((found nil))
            (while (and (not found) (not (bobp)))
              (forward-line -1)
              (let ((prev (gmic--current-line-string)))
                (unless (string-match-p "^\\s-*\\(?:#.*\\)?$" prev)
                  (if in-string
                      ;; Inside a math string: base indent on paren delta of
                      ;; the previous line.
                      ;; We add back paren-leading-close-delta of the previous
                      ;; line because those leading closers were already
                      ;; "spent" to dedent the previous line itself — they
                      ;; must not reduce the next line's indent a second time.
                      ;; Example: after "    );" (indent=4, delta=-2, leading=2)
                      ;; the next line should sit at 4+(-2)+2=4, not 4+(-2)=2.
                      (progn
                        (setq indent (+ (gmic--line-indentation prev)
                                        (gmic--paren-delta prev)
                                        (gmic--paren-leading-close-delta prev)))
                        (when (gmic--line-opens-string-p prev)
                          (setq indent (+ indent gmic-indent-offset)))
                        (setq found t))
                    ;; Outside strings: skip lines that are inside a string.
                    (let ((prev-in-string (gmic--in-multiline-string-p)))
                      (unless prev-in-string
                        (setq indent (+ (gmic--line-indentation prev)
                                        (gmic--line-delta prev)))
                        (setq found t))))))))))))
    ;; A command definition is always at column 0.
    ;; (Skipped when we already handled a lone closing-quote line above.)
    (unless (and in-string (gmic--line-is-closing-quote-p current-line))
      (when (or (gmic--line-is-command-def-p current-line)
                (gmic--line-is-cli-gui-comment-p current-line))
        (setq indent 0))
      ;; Adjust for leading closers on the current line.
      ;; The matching-open path already returns the final column, so skip the
      ;; subtraction in that case.
      (if in-string
          (unless (gmic--find-matching-open-indent current-line)
            (setq indent (- indent (gmic--paren-leading-close-delta current-line))))
        (setq indent (- indent (gmic--line-leading-close-delta current-line))))
      ;; Clamp to zero.
      (setq indent (max 0 indent))
      (indent-line-to indent))))

;;;; -------------------------------------------------------------------------
;;;; Imenu support — navigate to command definitions

(defvar gmic-imenu-generic-expression
  '(("Commands" "^\\s-*\\([a-zA-Z_][a-zA-Z0-9_]*\\)\\s-*\\(([^)]*)\\s-*\\)?:\\([^=]\\|$\\)" 1))
  "Imenu expressions for `gmic-mode'.")

;;;; -------------------------------------------------------------------------
;;;; Run support

(defun gmic-run ()
  "Run the current G'MIC script buffer using `gmic-executable'."
  (interactive)
  (let ((file (buffer-file-name)))
    (if file
        (compile (concat gmic-executable " " (shell-quote-argument file)))
      (user-error "Buffer has no associated file; save it first"))))

(defun gmic-run-region (beg end)
  "Run the G'MIC commands in the selected region BEG to END."
  (interactive "r")
  (let ((code (buffer-substring-no-properties beg end))
        (tmpfile (make-temp-file "gmic-region" nil ".gmic")))
    (with-temp-file tmpfile (insert code))
    (compile (concat gmic-executable " " (shell-quote-argument tmpfile)))))

;;;; -------------------------------------------------------------------------
;;;; Keymap

(defvar gmic-mode-map
  (let ((map (make-sparse-keymap)))
    (define-key map (kbd "C-c C-c") #'gmic-run)
    (define-key map (kbd "C-c C-r") #'gmic-run-region)
    (define-key map (kbd "C-c C-d") #'gmic-lookup-doc)
    map)
  "Keymap for `gmic-mode'.")

;;;; -------------------------------------------------------------------------
;;;; Documentation lookup

(defun gmic-lookup-doc ()
  "Open G'MIC online documentation for the command at point."
  (interactive)
  (let ((word (thing-at-point 'symbol t)))
    (if word
        (browse-url (concat "https://gmic.eu/reference/" word ".html"))
      (browse-url "https://gmic.eu/reference/"))))

;;;; -------------------------------------------------------------------------
;;;; Mode definition

;;;###autoload
(define-derived-mode gmic-mode prog-mode "G'MIC"
  "Major mode for editing G'MIC script files.

Key bindings:
\\{gmic-mode-map}"
  :syntax-table gmic-mode-syntax-table

  ;; Comments
  (setq-local comment-start "# ")
  (setq-local comment-end "")
  (setq-local comment-start-skip "#+\\s-*")

  ;; Font-lock
  (setq-local font-lock-defaults
              '(gmic-font-lock-keywords
                nil   ; no keywords-only (so strings are highlighted)
                nil   ; case sensitive
                nil   ; no syntax-alist override
                nil)) ; no syntax-begin override

  ;; Indentation
  (setq-local indent-line-function #'gmic-indent-line)

  ;; Imenu
  (setq-local imenu-generic-expression gmic-imenu-generic-expression)
  (imenu-add-to-menubar "Commands")

  ;; Misc
  (setq-local require-final-newline t)
  (setq-local indent-tabs-mode nil)
  (setq-local tab-width gmic-indent-offset))

;;;###autoload
(add-to-list 'auto-mode-alist '("\\.gmic\\'" . gmic-mode))

(provide 'gmic-mode)
;;; gmic-mode.el ends here
