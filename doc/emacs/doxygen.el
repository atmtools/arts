;;; doxygen.el --- support for doxygen style comments

;; Copyright (C) 2000 Basis Technology, Corp.

;; Author: Tom Emerson <tree@basistech.com>
;; Created: 2000-07-04
;; Keywords: languages comments doxygen

;;; Commentary:

;; TODO:
;;
;; - better documentation
;; - key bindings, perhaps
;; - allow field insertion, a la BibTeX mode
;; - generalize comment types - right now these reflect my personal
;;   style and the fact that I'm doing all of my work in C++.

;; Stefan Buehler 2000-09-13
;; - Replace *! by **. This means that always the first sentence is
;;   used as the short description.  
;; - Removed dummy text, takes only time to remove this.

;;; Code:

(defvar doxygen-date-format "%Y-%m-%d"
  "The format used to display dates when using the \\date command.")

(defun doxygen-insert-comment()
  "Insert a generic Doxygen comment block at point, including brief
and long sections."
  (interactive "*")
  (insert (concat "/**\n"
                  "   \n"
                  "*/\n"))
  (re-search-backward " ")
  (end-of-line))

(defun doxygen-insert-file-comment ()
  "Insert a Doxygen file comment at point."
  (interactive "*")
  (let ((file-name (if (buffer-file-name)
                       (file-name-nondirectory (buffer-file-name))
                     "untitled"))
        (date-string (format-time-string doxygen-date-format))
        (who (user-full-name)))
    (insert (format (concat "/**\n"
                            "   \\file   %s\n"
                            "\n"
                            "   \n"
                            "\n"
                            "   \\author %s\n"
                            "   \\date   %s\n"
                            "*/\n")
                    file-name who date-string)))
  (re-search-backward "   ")
  (re-search-backward "   ")
  (re-search-backward "   ")
  (re-search-backward "   ")
  (end-of-line))


(defun doxygen-insert-function-comment ()
  "Insert a Doxygen comment for the function at point."
  (interactive "*")
  (beginning-of-line)
  (let ((args (find-arg-list))
	(date-string (format-time-string doxygen-date-format))
        (who (user-full-name)))
    (insert (concat "/**\n"
                    "   \n"
                    "\n"))
    (when (cdr (assoc 'args args))
         (dump-arguments (cdr (assoc 'args args))))
    (unless (string= "void" (cdr (assoc 'return args)))
      (insert "   \\return \n"))
    (insert (format (concat "\n"
                            "   \\author %s\n"
                            "   \\date   %s\n"
                            "*/\n")
                    who date-string)))
  (re-search-backward "/")
  (re-search-backward "/")
  (re-search-forward " ")
  (end-of-line))

(defun doxygen-insert-member-group-region (start end)
  "Make the current region a member group."
  (interactive "*r")
  (save-excursion
    (goto-char start)
    (beginning-of-line)
    ; indent-according-to-mode doesn't work well here...
    (insert "//@{\n")
    (goto-char end)
    (end-of-line)
    (insert "\n//@}\n")))

(defun doxygen-insert-compound-comment ()
  "Insert a compound comment."
  (interactive "*")
  (let ((comment-start "//!< ")
        (comment-end ""))
    (indent-for-comment)))


;;; internal utility functions

(defun dump-arguments (arglist)
  "Insert a comment with the Doxygen comments for a function."
  (mapcar (function (lambda (x)
                      (insert (format "   \\param %s\n"
                                      (extract-argument-name x)))))
          arglist))

(defun extract-argument-name (arg)
  "Get the argument name from the argument string 'arg'."
  ; this does *not* work for function pointer arguments
  (if (string-match "\\(\\w+\\)\\(\\[\\]\\)?$" arg)
      (substring arg (match-beginning 1) (match-end 1))
    arg))

(defun find-arg-list ()
  "Extract various bits of information from a C or C++ function declaration"
  (save-excursion
    (if (re-search-forward (concat
                            ; want to use a shy group here, but I cannot
                            ; get it to work: "\\(?:const \|static \\)*"
                            ; fails to do the Right Thing. Urgh.
                            "\\(virtual \|const \|static \\)*?" ; opt. modifiers
                            "\\([a-zA-Z0-9_]+\\) "    ; return type
                            "\\([a-zA-Z0-9_:]+\\)"    ; function name
                            "(\\([^)]*\\))"           ; argument list
                            "\\([ \t]+const\\)?")     ; possibly const
                           nil t)
        (list (cons 'return   (buffer-substring (match-beginning 2)
                                                (match-end 2)))
              (cons 'function (buffer-substring (match-beginning 3)
                                                (match-end 3)))
              (cons 'args     (split-string
                               (buffer-substring (match-beginning 4)
                                                 (match-end 4)) ", +")))
      nil)))

(provide 'doxygen)

;;; doxygen.el ends here

