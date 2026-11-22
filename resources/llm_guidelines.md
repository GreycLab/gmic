# A Reminder of Best Practices for Coding G'MIC Filters

## General rules
1. Use English for the names of filters, parameters, variables, and comments in the code.
2. Use two spaces for each level of indentation.
3. If you need to test the execution of your code with a recent version of G'MIC, install the G'MIC binding for Python (`pip install gmic`).
4. No line must exceed 120 characters, indentation included. Use `\` to break long-lines into multiple lines.

## Command pipeline
5. In a G'MIC pipeline, a blank character acts as a separator between items. Never include spaces in list of arguments or variable assignments. E.g., `var = 3` is incorrect, but `var=3` is correct. If spaces are needed in argument (e.g. a string), embrace the argument with double-quotes.
6. To process each image in the list individually, use the `foreach { ... }` construct instead of `repeat $! l[$>] ... done done`.
7. To apply a command to an existing image and get the result a new image at the end of the list, use the prefix `+` before the command name. E.g., `+mirror[0] x`. If you need to clone an existing image and place it at the end of the list, simply use `[n]`, where `n` is the index or name of the image to be cloned (command `dup` does not exist).
8. Do not use `{}` around an expression evaluated as part of a control statement that evaluates to a Boolean. E.g., write `if round(x)>=pi` rather than `if {round(x)>=pi}`.
9. When assigning a value to a variable that requires a calculation, use `var:=` rather than `var={...}`. E.g., `x:=round(pi)` is preferable.
10. In recent G'MIC versions, command `endl` no longer exists; it has been replaced by `done` (or `}`, which is equivalent).
11. Use command shortcuts when possible (e.g. `r` for `resize`, `f` for `fill`, `c` for `cut`, `n` for `normalize`, `s` for `split`, etc.)
12. Use selection shortcut when possible (e.g. `.` for `[-1]`, `..` for `[-2]` and `...` for `[-3]`).
13. To create a new image filled with a mathematical expression `"expr"`, use `width,height,depth,spectrum,"expr"` instead of `width,height,depth,spectrum f. "expr"`.

## Mathematical expressions
14. When you specify mathematical expressions enclosed in double quotes (e.g. while using the `fill` and `eval` commands), you are allowed to write expressions across multiple lines and include spaces between operators. Use these features to make your code clearer and more readable.
15. In a mathematical expression, if you need to divide a term by `100` , use the `%` suffix instead of writing `/100`. E.g., write `x%` instead of `x/100`.
16. When assigning a variable in a math expression, put spaces around the `=` sign, e.g. `x = 0` rather than `x=0`.
17. When adding or subtracting two terms, put spaces around the `+` or `-` operator, e.g. `x + y`, rather than `x+y`.
18. Use `°` rather than `*pi/180`, e.g. `x = angle°` rather than `x = angle*pi/180`.

## Definition of filter parameters and GUI
19. When defining a G'MIC-Qt filter, the main command name must begin with the prefix `fx_`. If a command to generate the preview is required, it must have the same name, suffixed by `_preview`.
20. When defining the filter parameter GUI, the three first lines after the definition must be:
~~~
#@gui : Description = note{"<span color="#EE5500"><b>Description:</b></span> \
# A one-line description of what the filter does. Keep it short and concise.
#@gui : "}, sep = separator()
~~~
21. When defining the filter parameter GUI, the two last lines must be:
~~~
#@gui : sep = separator()
#@gui : note = note("<small>Author: <i>Author Name</i>.      Latest Update: <i>YYYY/MM/DD</i>.</small>")
~~~
where `YYYY/MM/DD` is the date of the day. There must be exactly 6 non-breaking spaces before `Latest Update:`.
22. Parameter names must be capitalized as in an English title, following APA conventions
23.  Do not include a `Random Seed` parameter for the filter unless you are certain that the rendering displayed in the preview window will be roughly the same as the final rendering once the filter is applied.
24. When defining a `choice()` parameter, if the default choice is the first item, you don't need to specify a `0` as the first argument.
25. If it might be useful to randomize a parameter so that the user can randomly explore the filter's parameter space, prefix the parameter type with `~`. E.g., `Color = ~color(255,255,255)`.
26. If certain filter parameters share a common theme (such as "Colors" or "Geometry"), group them between two `separator()`, and add a group header note after the first separator, such as:
~~~
#@gui : note = note("<small><span color="#CC3300"><b>Colors:</b></span></small>")
~~~
27. If a filter has more than 6 parameters, define variables for each of these parameters, e.g. with:
~~~
seed,randomize,nbd,dsz,imode,istr,flo,swi,det,palsel,r1,g1,b1,r2,g2,b2,pr,pg,pb=${1-19}
~~~
then use these named variables in the code, rather than the function’s direct arguments (`$1`, `$2`, etc.), for better readability.
28. In the code for a G'MIC-Qt filter, do not specify default argument values (e.g., `${1=10},${2=20}`). All arguments are guaranteed to get a value when the filter is called from the plug-in.
