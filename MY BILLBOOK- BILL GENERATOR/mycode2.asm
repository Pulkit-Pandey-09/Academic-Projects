.model small
.stack 100h
.data
menu db '*****************FOOD MENU*********************$'
menu1 db 'Press 1 for TEA$ '
menu2 db 'Press 2 for COFFEE$'
menu3 db 'Press 3 for BURGER$'  
menu4 db 'Press 4 to COCA COLA $'
menu6 db 'Press 6 to FRIES$'  
menu5 db 'Press 5 to TACOS$'
menu7 db 'Press 7 to SHOW THE RECORD$'
menu8 db 'Press 8 to DELETE THE RECORD$'
menu9 db 'Press 9 to EXIT$' 

msg1 db ' OUT OF STOCK$'
msg2 db 'WRONG INPUT$'
msg5 db 'RECORD$'
msg6 db 'YOU CAN ORDER MORE$'
msg7 db 'THE TOTAL AMOUNT IS=$'
msg8 db 'THE TOTAL NUMBER OF ITEMS ORDERED=$'
msg9 db 'TOTAL NUMBER OF TEA ORDERED=$'
msg10 db 'TOTAL NUMBER OF COFFEE ORDERED=$'
msg11 db 'TOTAL NUMBER OF BURGER ORDERED=$'
msg12 db '***      RECORD DELETED SUCCESSFULLY       ***$'  
msg13 db 'TOTAL NUMBER OF COCO COLA ORDERED=$' 
msg14 db 'TOTAL NUMBER OF FRIES ORDERED=$'
msg15 db 'TOTAL NUMBER OF TACOS ORDERED=$'
amount dw 0                                
count dw  '0'
am1 dw ?
am2 dw ?
am3 dw ?


r db '0'
c db '0'
b db '0' 
x db '0' 
y db '0'
z db '0'

.code
main proc
mov ax,@data
mov ds,ax

;mov cx,count

;mov cx,0

while_:   
             ;Menu
mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h


mov dx,offset menu
mov ah,9
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h


mov dx,offset menu1
mov ah,9
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h


mov dx,offset menu2
mov ah,9
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h


mov dx,offset menu3
mov ah,9
int 21h
mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h


mov dx,offset menu4
mov ah,9
int 21h
mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h



mov dx,offset menu5
mov ah,9
int 21h
mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h


mov dx,offset menu6
mov ah,9
int 21h
mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h

mov dx,offset menu7
mov ah,9
int 21h
mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h

mov dx,offset menu8
mov ah,9
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h

mov dx,offset menu9
mov ah,9
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h


  ;userinput

mov ah,1
int 21h
mov bl,al
mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h

  ;now compare
mov al,bl
cmp al,'1'
je rikshw
cmp al,'2'
je car
cmp al,'3'
je bus
cmp al,'7'
je rec
cmp al,'8'
je del
cmp al,'9'
je end_  
cmp al,'6'
je xyz
cmp al,'5'
je xyz1
cmp al,'4'
je xyz2

mov dx,offset msg2
mov ah,9
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h
jmp while_

rikshw:
call rikshaw


car:
call caar


rec:
call recrd


del:
call delt


bus:
call buss

xyz:
call xyzx   

xyz1:
call xyzy

xyz2:
call xyzz



end_:
mov ah,4ch
int 21h

main endp





rikshaw proc
cmp count,'8'
jle rikshw1
mov dx,offset msg1
mov ah,9
int 21h
jmp while_
jmp end_

rikshw1:
mov ax,20
add amount, ax
mov dx,0 ; remainder is 0
mov bx,10 
mov cx,0
l2:
        div bx
        push dx
        mov dx,0
        mov ah,0
        inc cx
        cmp ax,0
        jne l2
   
l3:
        pop dx
        add dx,48
        mov ah,2
        int 21h
loop l3
;mov am1,dx
inc count
;mov dx,count
inc r

jmp while_
jmp end_

xyzx proc
cmp count,'8'
jle xyzxx
mov dx,offset msg1
mov ah,9
int 21h
jmp while_
jmp end_

xyzxx:
mov ax,50
add amount, ax
mov dx,0 ; remainder is 0
mov bx,10 
mov cx,0
lxyzx2:
        div bx
        push dx
        mov dx,0
        mov ah,0
        inc cx
        cmp ax,0
        jne lxyzx2
   
lxyzx3:
        pop dx
        add dx,48
        mov ah,2
        int 21h
loop lxyzx3
;mov am1,dx
inc count
;mov dx,count
inc x

jmp while_
jmp end_

xyzy proc
cmp count,'8'
jle xyzyy
mov dx,offset msg1
mov ah,9
int 21h
jmp while_
jmp end_

xyzyy:
mov ax,40
add amount, ax
mov dx,0 ; remainder is 0
mov bx,10 
mov cx,0
lxyzy2:
        div bx
        push dx
        mov dx,0
        mov ah,0
        inc cx
        cmp ax,0
        jne lxyzy2
   
lxyzy3:
        pop dx
        add dx,48
        mov ah,2
        int 21h
loop lxyzy3
;mov am1,dx
inc count
;mov dx,count
inc y

jmp while_
jmp end_ 
     

xyzz proc
cmp count,'8'
jle xyzzz
mov dx,offset msg1
mov ah,9
int 21h
jmp while_
jmp end_

xyzzz:
mov ax,35
add amount, ax
mov dx,0 ; remainder is 0
mov bx,10 
mov cx,0
lxyzz2:
        div bx
        push dx
        mov dx,0
        mov ah,0
        inc cx
        cmp ax,0
        jne lxyzz2
   
lxyzz3:
        pop dx
        add dx,48
        mov ah,2
        int 21h
loop lxyzz3
;mov am1,dx
inc count
;mov dx,count
inc z

jmp while_
jmp end_ 
 
     
caar proc
cmp count,'8'
jle car1
mov dx,offset msg1
mov ah,9
int 21h
jmp while_
jmp end_

car1:
mov ax,30
add amount, ax
mov dx,0
mov bx,10
mov cx,0
l22:
        div bx
        push dx
        mov dx,0
        mov ah,0
        inc cx
        cmp ax,0
       jne l22
   
l33:
        pop dx
        add dx,48
        mov ah,2
        int 21h
loop l33

;mov am2,amount

inc count
inc c
jmp while_
jmp end_




buss proc
cmp count,'8'
jle bus1
mov dx,offset msg1
mov ah,9
int 21h
jmp while_
jmp end_

bus1:
mov ax,50
add amount, ax
mov dx,0
mov bx,10
mov cx,0
l222:
        div bx
        push dx
        mov dx,0
        mov ah,0
        inc cx
        cmp ax,0
       jne l222
   
l333:
        pop dx
        add dx,48
        mov ah,2
        int 21h
loop l333
;mov am3,amount

inc count
inc b
jmp while_
jmp end_


recrd proc
mov dx,offset msg7
mov ah,9
int 21h


; print here the whole amount
mov ax, amount

mov dx,0
mov bx,10
mov cx,0
totalpush:
        div bx
        push dx
        mov dx,0
      ;  mov ah,0
        inc cx
        cmp ax,0
       jne totalpush
   
totalprint:
        pop dx
        add dx,48
        mov ah,2
        int 21h
loop totalprint




mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h







mov dx,offset msg8
mov ah,9
int 21h

mov dx,count
mov ah,2
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h


mov dx,offset msg9
mov ah,9
int 21h

mov dl,r
mov ah,2
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h

mov dx,offset msg13
mov ah,9
int 21h

mov dl,z
mov ah,2
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h
       

mov dx,offset msg14
mov ah,9
int 21h

mov dl,x
mov ah,2
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h


mov dx,offset msg15
mov ah,9
int 21h

mov dl,y
mov ah,2
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h
       

mov dx,offset msg10
mov ah,9
int 21h


mov dl,c
mov ah,2
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h



mov dx,offset msg11
mov ah,9
int 21h

mov dl,b
mov ah,2
int 21h

jmp while_
jmp end_


delt proc
mov r,'0'
mov c,'0'
mov b,'0' 
mov x,'0'
mov y,'0'
mov z,'0'

mov amount,0
;sub amount,48
mov count,'0'
mov dx,offset msg12
mov ah,9
int 21h

mov dx,10
mov ah,2
int 21h
mov dx,13
mov ah,2
int 21h




jmp while_
jmp end_