FUNCTION colored, str, c=c, s=s
    ;+#colored
    ;+Creates colored string
    ;+***
    ;+##Arguments
    ;+     **str**: String to be colored
    ;+
    ;+##Keyword Arguments
    ;+     **c**: Foreground color code
    ;+
    ;+     **s**: Style code
    ;+
    ;+###Foreground Color Codes
    ;+     **k**: Black,
    ;+     **r**: Red,
    ;+     **g**: Green,
    ;+     **y**: Yellow,
    ;+     **b**: Blue,
    ;+     **m**: Magenta,
    ;+     **c**: Cyan,
    ;+     **w**: White
    ;+###Style Format Codes
    ;+     **n**: Normal, 
    ;+     **b**: Bright, 
    ;+     **d**: Dim, 
    ;+     **i**: Italics, 
    ;+     **u**: Underline, 
    ;+     **r**: Reverse, 
    ;+     **h**: Hidden, 
    ;+     **s**: Strikethrough 
    ;+##Return Value
    ;+     Colored string
    ;+
    ;+##Example Usage
    ;+```idl
    ;+IDL> blue_bright_hello = colored("Hello",c="b",s="b")
    ;+```
    if not keyword_set(c) then c='w' ;;Foreground Color
    if not keyword_set(s) then s='n' ;;Style
    esc=string(27b)
    back=esc+"[0m"

    style = {n:'0',b:'1',d:'2',i:'3',u:'4',r:'7',h:'8',s:'9'}
    sTags = ["n","b","d","i","u","r","h","s"]

    fgColors = { k:'30',r:'31',g:'32',y:'33',$
                b:'34',m:'35',c:'36',w:'37'}
    fgTags = ["k","r","g","y","b","m","c","w"]

    sIndex = where(s eq sTags)
    fgIndex = where(c eq fgTags)
    if sIndex eq -1 then sIndex=0
    if fgIndex eq -1 then fgIndex=7

    return, esc+"["+style.(sIndex)+";"+fgColors.(fgIndex)+"m"+str+back
END
