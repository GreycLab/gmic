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
  ;; Commandes natives G'MIC, hard-codées en C++ — liste stable.
  ;; Les opérateurs purement symboliques (!=, %, &, *, +, -, /, etc.)
  ;; sont gérés séparément via font-lock-operator-face.
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

(defvar gmic--math-functions
  '("abs" "acos" "arg" "argkth" "argmax" "argmin" "argmaxabs" "argminabs"
    "asin" "atan" "atan2" "avg" "begin" "begin_t" "bool" "cbrt" "ceil"
    "complex" "conj" "copy" "correlate" "cos" "cosh" "critical" "cross"
    "cut" "da_back" "da_capacity" "da_find" "da_freeze" "da_insert"
    "da_pop" "da_push" "da_remove" "da_resize" "da_size" "date" "debug"
    "det" "diag" "dot" "e" "eig" "ellipse" "end" "end_t" "eval" "exp"
    "expr" "eye" "f2ui" "fact" "find" "floor" "fsize" "gauss"
    "haar" "ihaar" "id" "if" "ilog2" "im" "inf" "inrange" "int"
    "inv" "isfile" "isdir" "isinf" "isnan" "isnum" "isvec"
    "kth" "lerp" "linear" "log" "log2" "log10" "logit" "lowercase"
    "lu" "map" "max" "maxabs" "median" "merge" "min" "minabs"
    "mul" "nan" "norm" "normP" "norminf" "normalize" "o" "permute"
    "pi" "print" "prod" "pseudoinv" "q" "qr" "rand" "re" "resize"
    "reverse" "rol" "ror" "rot" "round" "run"
    "s" "same" "set" "sign" "sin" "sinc" "sinh" "size" "solve" "sort"
    "sqr" "sqrt" "std" "stod" "store" "strcat" "strcmp" "strcopy"
    "strpbrk" "strrstr" "strstr" "strtod" "strtol" "strupr"
    "sum" "svd" "tan" "tanh" "trace" "transpose" "trunc" "type" "ui2f"
    "unitnorm" "uppercase" "var" "vector" "vmax" "vmin"
    "whiledo" "xor")
  "G'MIC math expression functions.")

(defconst gmic-font-lock-keywords
  (list
   ;; 1. Command definitions: "name :" or "name(args) :"
   ;;    Le ':' ne doit pas être suivi de '=' (sinon c'est une affectation var:=)
   ;;    On utilise :[^=] ou : en fin de ligne (pas de lookahead en Emacs regexp)
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

   ;; 7. Math functions inside double-quoted strings or braces
   (cons (concat "\\_<"
                 (regexp-opt gmic--math-functions)
                 "\\s-*(")
         '(0 font-lock-type-face))

   ;; 8. Native built-in commands — anywhere on the line, at word boundary
   (cons (concat "\\_<" (regexp-opt gmic--builtin-commands t) "\\_>")
         'font-lock-builtin-face)

   ;; 9. Symbolic native operators: +3d -3d *3d /3d m* m/ => != == <= >= << >>
   '("\\(?:\\+3d\\|-3d\\|\\*3d\\|/3d\\|m\\*\\|m/\\|!=\\|==\\|<=\\|>=\\|<<\\|>>\\|=>\\)"
     . font-lock-builtin-face)

   ;; 10. Numeric literals (integers and floats, including degree notation)
   '("\\b\\([0-9]+\\.?[0-9]*\\(?:e[+-]?[0-9]+\\)?°?\\)\\b"
     (1 font-lock-constant-face))

   ;; 11. Shebang line
   '("^#!.*$" . font-lock-comment-face)

   ;; 12. Special constants
   '("\\_<\\(pi\\|inf\\|nan\\|true\\|false\\)\\_>" . font-lock-constant-face))
  "Font-lock keywords for `gmic-mode'.")

;;;; -------------------------------------------------------------------------
;;;; Indentation — comptage des ouvrants/fermants sur chaque ligne

(defvar gmic--indent-open-re
  ;; Seuls les mots-clés comptent comme ouvrants.
  ;; '{' n'est jamais compté : c'est un raccourci syntaxique pour le mot-clé
  ;; qui le précède, lequel est déjà compté.
  (concat "\\_<"
          (regexp-opt '("repeat" "for" "foreach" "do"
                        "if" "elif" "else"
                        "local" "l")
                      t)
          "\\_>")
  "Regexp matching block-opening keywords in G'MIC.")

(defvar gmic--indent-close-re
  ;; '}' isolé (précédé d'un blanc ou en début de ligne, suivi d'un blanc ou fin)
  ;; + mots-clés fermants.
  (concat "\\(?:\\(?:^\\|\\s-\\)}\\(?:\\s-\\|$\\)\\|\\_<"
          (regexp-opt '("done" "fi" "elif" "else" "while") t)
          "\\_>\\)")
  "Regexp matching block-closing keywords or isolated '}' in G'MIC.")

(defun gmic--strip-comments-and-strings (line-str)
  "Return LINE-STR with comments and string contents removed.
Keeps the quote delimiters themselves to preserve quote counting,
but removes everything between them."
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
         ((and (not in-string) (= ch ?#))
          (setq i len))               ; stop — rest is comment
         ((not in-string)
          (setq result (concat result (string ch))))))
      (setq i (1+ i)))
    result))

(defun gmic--strip-comments (line-str)
  "Return LINE-STR with any trailing G'MIC comment removed.
A comment starts at '#' unless it is inside a double-quoted string.
String contents are preserved (unlike `gmic--strip-comments-and-strings')."
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
         ((and (not in-string) (= ch ?#))
          (setq i len))
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
         ;; Fermants sur la même ligne qu'un ouvrant annulent cet ouvrant.
         ;; Ex: "repeat 6 { ... done" sur une ligne = delta 0.
         ;; On ne compte que les fermants qui suivent un ouvrant sur la ligne.
         (close-re (concat "\\(?:\\(?:^\\|\\s-\\)}\\(?:\\s-\\|$\\)\\|\\_<"
                           (regexp-opt '("done" "fi" "while") t)
                           "\\_>\\)"))
         ;; Fermants qui apparaissent APRÈS le premier ouvrant de la ligne
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
         (leading-elif-else
          (if (string-match-p
               (concat "^\\s-*\\_<\\(?:elif\\|else\\)\\_>")
               stripped)
              1 0))
         ;; Autres fermants : done, fi, while, '}' isolé — avant le premier ouvrant
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
          (setq i len)))) ; fin de ligne : commentaire
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
      ;; Remonter jusqu'à la définition de commande englobante ou le début du buffer
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

(defun gmic--line-opens-string-p (line-str)
  "Return t if LINE-STR opens a multiline string (odd number of quotes)."
  (= (% (gmic--count-unescaped-quotes line-str) 2) 1))

(defun gmic-indent-line ()
  "Indent the current line for `gmic-mode'.

Algorithm:
  1. Determine if the current line is inside a multiline string.
  2. Find the previous non-empty line that is NOT inside a multiline
     string (i.e. the last structural reference line).
  3. Start from that line's indentation + delta.
  4. If we are inside a multiline string, add one extra level.
  5. Subtract the leading-close contribution of the CURRENT line.
  6. Clamp to zero."
  (interactive)
  (let* ((in-string (gmic--in-multiline-string-p))
         (indent 0))
    (save-excursion
      (beginning-of-line)
      (let ((found nil))
        (while (and (not found) (not (bobp)))
          (forward-line -1)
          (let ((prev (gmic--current-line-string)))
            (unless (string-match-p "^\\s-*\\(?:#.*\\)?$" prev)
              ;; N'utiliser cette ligne comme référence que si elle
              ;; n'est pas elle-même à l'intérieur d'une chaîne
              (let ((prev-in-string (gmic--in-multiline-string-p)))
                (unless prev-in-string
                  (setq indent (+ (gmic--line-indentation prev)
                                  (gmic--line-delta prev)))
                  (setq found t))))))))
    ;; Indentation supplémentaire si on est dans une chaîne multiligne
    (when in-string
      (setq indent (+ indent gmic-indent-offset)))
    ;; Une définition de commande est toujours à la colonne 0
    (when (gmic--line-is-command-def-p (gmic--current-line-string))
      (setq indent 0))
    ;; Adjust for leading closers on the current line
    (let ((cur (gmic--current-line-string)))
      (setq indent (- indent (gmic--line-leading-close-delta cur))))
    ;; Clamp
    (setq indent (max 0 indent))
    (indent-line-to indent)))

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
