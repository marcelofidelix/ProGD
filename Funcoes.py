from numpy import radians, cos, sin, argmax, round_

 #Função que permite o cálculo do raio de ação
def raio(teta, L, J, S):
	teta_rad = radians(teta)
	return  J + L*cos(teta_rad) + S*sin(teta_rad)

#CÁLCULO DO PARÂMETRO Vd. O PARÂMETRO Vd CORRESPONDE À VELOCIDADE DO
#BARCO. OS VALORES SÃO CALCULADOS CONFORME TABELA 3, EM ft/s
def calc_Vd(esc_embarc,Hsig):
    Vd = 0
    if (esc_embarc == 'Estrutura fixa'):
        Vd = 0. #ft/s
    elif (esc_embarc == 'Embarcação' and Hsig < 9.8):
        Vd = .6*Hsig #ft/s
    elif (esc_embarc == 'Embarcação' and Hsig >= 9.8):
        Vd = 5.9 + .3*(Hsig - 9.8) #ft/s
    return Vd
 
#CÁLCULO DO PARÂMETRO Vc. O PARÂMETRO Vc CORRESPONDE À VELOCIDADE DA
#PONTA DA LANÇA DO GUINDASTE, DEVIDO À MOVIMENTAÇÃO DA UEP. É OBTIDO
#ATRAVÉS DA TABELA 3 EM ft/s
def calc_Vc(esc_UEP,Hsig):
    Vc = 0
    if (esc_UEP == 'Fixa'):
        Vc = 0. #ft/s
    elif (esc_UEP == 'SS'):
        Vc = 0.025 * Hsig * Hsig
    elif (esc_UEP == 'FPSO'):
        Vc = 0.05 * Hsig * Hsig
    return Vc

#CÁLCULO DO PARÂMETRO Av. O PARÂMETRO Av CORRESPONDE À ACELERAÇÃO
#VERTICAL DO GUINDASTE. A UNIDADE É *g, ONDE g = 32,2ft/s²

def calc_Av(esc_UEP,Hsig):
    if (esc_UEP == 'Fixa'):
        Av = 0. #*g
    elif (esc_UEP == 'SS'):
        Av = 0.0007 * Hsig * Hsig
        if (Av < .07):
            Av = .07
    elif (esc_UEP == 'FPSO'):
        Av = 0.0012 * Hsig * Hsig
        if (Av < .07):
            Av = .07
    return Av

#CHA CORRESPONDE AO CRANE DYNAMIC HORIZONTAL ACCELERATION, CALCULADO
#CONFORME A TABELA 5. UNIDADE: *g

def calc_CHA(esc_UEP,Hsig):
    CHA = ''
    if (esc_UEP == 'Fixa'):
        CHA = 0. #*g
    elif (esc_UEP == 'SS'):
        CHA = 0.007 * Hsig
        if (CHA < .03):
            CHA = .03
    elif (esc_UEP == 'FPSO'):
        CHA = .01*(Hsig)**1.1
        if (CHA < .03):
            CHA = .03
    return CHA

#OS PARÂMETROS LIST E TRIM CORRESPONDEM AO ÂNGULO DE INCLINAÇÃO DA
#BASE DOS GUINDASTES. VALORES EM °

def calc_List(esc_UEP):
    List = 0
    if (esc_UEP == 'Fixa'):
        List = .5 #°
    elif (esc_UEP == 'SS'):
        List = 1.5
    elif (esc_UEP == 'FPSO'):
        List = 2.5
    return List

def calc_Trim(esc_UEP):
    Trim = 0
    if (esc_UEP == 'Fixa'):
        Trim = .5 #°
    elif (esc_UEP == 'SS'):
        Trim = 1.5
    elif (esc_UEP == 'FPSO'):
        Trim = 1.
    return Trim

#PARÂMETRO Vhmin: During offboard lifts, hoisting velocity at the
#elevation where the lift is initiated (i.e. supply boat deck level)
#shall be fast enough to avoid re-contact after the load is lifted.
#The minimum hoisting velocity (Vhmin) for any particular hook load to
#be lifted shall be

def calc_Vhmin(Hsig):
    if (Hsig <= 6.):
        Vhmin = 0.033 + 0.098 * Hsig
    elif (Hsig > 6.):
        Vhmin = 0.067 * (Hsig + 3.3)
    return Vhmin
