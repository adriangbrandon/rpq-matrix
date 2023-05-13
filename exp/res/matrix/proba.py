import re

def eliminar_parentesis_innecesarios(expresion):
    # Eliminar paréntesis innecesarios
    while True:
        nueva_expresion = re.sub(r'\([^()]+\)', '', expresion)
        if nueva_expresion == expresion:
            break
        expresion = nueva_expresion
    
    # Eliminar paréntesis redundantes
    nueva_expresion = re.sub(r'\(([^()]+)\)', r'\1', expresion)
    
    return nueva_expresion

# Ejemplo de uso
expresion = "((3 + 5) * 2) - (4 / (2 + 1))"
resultado = eliminar_parentesis_innecesarios(expresion)
print(resultado)

