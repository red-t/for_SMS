### Germ train & test ###
def func(n):
    a = [str(n)]
    b = [str(n)]
    for i in range(2,4):
        a.append(str(int(n / i)))
        b.append(str(n))

    print(" ".join(a))
    print(" ".join(b))

### Soma train ###
def func(n):
    a = [str(n)]
    b = [str(n)]
    for i in range(5, 31, 5):
        a.append(str(int(n / i)))
        b.append(str(n))
    
    print(" ".join(a))
    print(" ".join(b))

### Soma test1 ###
def func(n):
    a = [str(n)]
    b = [str(n)]
    for i in range(5, 11, 5):
        a.append(str(int(n / i)))
        b.append(str(n))
    
    print(" ".join(a))
    print(" ".join(b))

### Soma test2 ###
def func(n):
    a = []
    b = []
    for i in [15, 20, 25, 30, 40, 50]:
        a.append(str(n))
        b.append(str(int(n * i)))
    
    print(" ".join(a))
    print(" ".join(b))