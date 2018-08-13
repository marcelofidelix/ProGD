import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input
import plotly.graph_objs as go
import pandas as pd
import ast
from numpy import radians, cos, array, arctan, arcsin, sin, pi, rad2deg, degrees, copy, clip
import numpy as np
import Funcoes

app = dash.Dash()

#Lê o banco de dados e gera um DataFrame
df_dados = pd.read_csv('df_dados.csv')

#Ordena os dados por ordem alfabética
df_dados.sort_values('Modelo',inplace=True)

#Cria um array vazio
combo_modelos_valores = []

#Popula o array com dicionários que irão aparecer no dropdown dos modelos
for i in df_dados['Modelo']:
    combo_modelos_valores.append({'label': i, 'value': i})

#Lista de graficos que irão habitar o dropdown com os parâmetros para plotagem
lista_param = ['[Gráfico] Ângulo alfa',
'[Gráfico] Ângulo beta',
'[Gráfico] Capacidade Estática',
'[Gráfico] Capacidade Estática Orig. x Otim.',
'[Gráfico] Raio',
'[Gráfico] Cabo do sistema principal',
'[Gráfico] Cabo do sistema principal Orig. x Otim.',
'[Gráfico] Cabo de lança',
'[Gráfico] Cabo de lança Orig. x Otim.',
'[Gráfico] cvon',
'[Gráfico] cvonot',
'[Gráfico] Esforço de compressão na lança',
'[Gráfico] Esforço de compressão na lança Orig. x Otim.',
'[Gráfico] Esforço nas hastes traseiras do cav. Orig. x Otim.',
'[Gráfico] Esforço nas hastes traseiras do cav.',
'[Gráfico] Esforço nas hastes dianteiras do cav. Orig. x Otim.',
'[Gráfico] Esforço nas hastes dianteiras do cav.',
'[Gráfico] Momento Orig. x Otim.',
'[Gráfico] Momento',
'[Gráfico] Sustentação da lança',
'[Gráfico] Sustentação da lança Orig. x Otim.',
'[Gráfico] alfacl',
'[Gráfico] Rp',
'[Gráfico] Rpx',
'[Gráfico] Rpy',         
'[Gráfico] gama_teta',
'[Gráfico] Tcgx',         
'[Gráfico] Tcgy',
'[Gráfico] Tclx',
'[Gráfico] Tcly',
'[Gráfico] gama',
'[Gráfico] teta',
'[Variável] a',
'[Variável] b',
'[Variável] V',
'[Variável] G',
'[Variável] L',
'[Variável] N',
'[Variável] S',
'[Variável] J',
'[Variável] Pbola',
'[Variável] Npl',
'[Variável] Nrl',
'[Variável] Npm',
'[Variável] Nrm',
'[Variável] Pplat',
'[Variável] Dplat',
'[Variável] Pcp',
'[Variável] Dcp',
'[Variável] tipocorda',
'[Variável] hl',
'[Variável] ll',
'[Variável] Dc',
'[Variável] t',
'[Variável] matcordas',
'[Variável] Tec',
'[Variável] Npend',
'[Variável] Dcabolan',
'[Variável] Dcabomoi',
'[Variável] Dcabopend',
'[Variável] Pl',
'[Variável] M',
'[Variável] Lpend',
'[Variável] Ljib',
'[Variável] Sjib',
'[Variável] CC1',
'[Variável] CC2',
'[Variável] CC3',
'[Variável] CC4',
'[Variável] D1',
'[Variável] D2',
'[Variável] D3',
'[Variável] D4',
'[Variável] Elan',
'[Variável] L',
'[Variável] Pmoi',
'[Variável] Vc',
'[Variável] CHA',
'[Variável] List',
'[Variável] Trim',
'[Variável] H',
'[Variável] Efm',
'[Variável] Efl',
'[Variável] FLFm',
'[Variável] FLFl',
'[Variável] tetac',
'[Variável] slider_penTcg'
]

sorted(lista_param)

#lista de dicionários que vai receber os valores que populam o dropdown
combo_param = []
for i in lista_param:
    combo_param.append({'label':i,'value':i})

#Lista de graficos que irão habitar o dropdown raixo x ângulo
combo_grafico = [
    {'label': 'Ângulo(°)', 'value': 'Ângulo(°)'},
    {'label': 'Raio(m)', 'value': 'Raio(m)'}
]

#Definição do layout do app
app.layout = html.Div([
    
    html.Div([
    html.H2('GDCalc'),

    html.H4('Modelo'),

    dcc.Dropdown(
    id='combo_modelos',
    options=combo_modelos_valores,
    value='American 5750R 90ft'
    ),

    ],style={'width': '45%', 'display': 'inline-block'}),

    html.Div([

    html.Div(
    
    [html.H4('Penalizações'),

    html.Div('Cabo do guincho'),

    dcc.Slider(
        id='slider_penTcg',
        min=0,
        max=100,
        step=1,
            marks={
        25: '25%',
        50: '50%',
        75 : '75%',
        100: '100%'
    },
        value=100,
    ),

    html.Div('Cabo da lança'),

    dcc.Slider(
        id='slider_penTcl',
        min=0,
        max=100,
        step=1,
            marks={
        25: '25%',
        50: '50%',
        75 : '75%',
        100: '100%'
    },
        value=100,
    ),

    html.Div('Comp. da Lança'),

    dcc.Slider(
        id='slider_penEcl',
        min=0,
        max=100,
        step=1,
            marks={
        25: '25%',
        50: '50%',
        75 : '75%',
        100: '100%'
    },
        value=100,
    ),

    html.Div('Momento'),

    dcc.Slider(
        id='slider_penMom',
        min=0,
        max=100,
        step=1,
            marks={
        25: '25%',
        50: '50%',
        75 : '75%',
        100: '100%'
    },
        value=100,
    ),

    html.Div('Hastes tras. do cav.'),

    dcc.Slider(
        id='slider_penFht',
        min=0,
        max=100,
        step=1,
            marks={
        25: '25%',
        50: '50%',
        75 : '75%',
        100: '100%'
    },
        value=100,
    ),

    html.Div('Hastes dian. do cav.'),

    dcc.Slider(
        id='slider_penFhd',
        min=0,
        max=100,
        step=1,
            marks={
        25: '25%',
        50: '50%',
        75 : '75%',
        100: '100%'
    },
        value=100,
    ),

    html.Div('Fator de penalização'),

    dcc.Slider(
        id='slider_penFator',
        min=0,
        max=100,
        step=1,
        marks={
        25: '25%',
        50: '50%',
        75 : '75%',
        100: '100%'
    },
        value=100,
    ),

    html.H4('Módulo de elasticidade do cabos'),

    html.Div('Moitão (N/mm²)'),

    dcc.Input(
        id='In_Ecm',
        value='99974',
        type='text',
    ),

    html.Div('Lança (N/mm²)'),

    dcc.Input(
        id='In_El',
        value='86874',
        type='text',
    ),

    html.Div('Pendentes (N/mm²)'),

    dcc.Input(
        id='In_Ep',
        value='86874',
        type='text',
    ),

    html.H4('Altura do chassi em relação ao mar (m)'),

    dcc.Slider(
        id='slider_H',
        min=0,
        max=50,
        step=1,
        value=20,
        marks={
        0: '0',
        10: '10',
        20: '20',
        30: '30',
        40: '40',
        50: '50'
        },
    ),

    html.H4('Altura de onda significativa'),

    dcc.Slider(
        id='slider_Hsig',
        min=0,
        max=5,
        step=0.1,
        value=2,
    ),

    html.Div(id='slider_Hsig_out'),

    html.H4('Vh (ft/s)'),

    dcc.Slider(
        id='slider_Vh',
        min=0,
        max=5,
        step=0.1,
        value=1.3,
    ),],style={'width': '80%'}),

    html.Div(id='slider_Vh_out'),

    html.H4('UEP'),

    dcc.RadioItems(
        id='radio_UEP',
        options=[
            {'label': 'Fixa', 'value': 'Fixa'},
            {'label': 'SS', 'value': 'SS'},
            {'label': 'FPSO', 'value': 'FPSO'},
        ],
        value='Fixa',
        labelStyle={'display': 'inline-block'}
    ),

    html.H4('Embarcação'),

    dcc.RadioItems(
        id='radio_embarc',
        options=[
            {'label': 'Embarcação', 'value': 'Embarcação'},
            {'label': 'Estrutura fixa', 'value': 'Estrutura fixa'},
        ],
        value='Embarcação',
        labelStyle={'display': 'inline-block'}
    ),

    html.H4('Mancal das roldanas'),

    dcc.RadioItems(
        id='radio_mancal',
        options=[
            {'label': 'Rolamento', 'value': 'Rolamento'},
            {'label': 'Bucha', 'value': 'Bucha'}
        ],
        value='Rolamento',
        labelStyle={'display': 'inline-block'}
    ),

    ],style={'width': '100%', 'display': 'inline-block','columnCount': 2}),
    
    html.H4('Parâmetros de saída'),

    dcc.Dropdown(
        id='combo_param',
        options=combo_param,
        value='[Gráfico] Capacidade Estática'
    ),

    dcc.RadioItems(
        id='combo_ang_raio',
        options=combo_grafico,
        value='Ângulo(°)',
        labelStyle={'display': 'inline-block'}
    ),

    html.Div(id='mostra_modelo')

])

@app.callback(
    Output('mostra_modelo', 'children'),
    [Input('combo_modelos', 'value'),
    Input('combo_ang_raio', 'value'),
    Input('combo_param', 'value'),
    Input('radio_mancal', 'value'),
    Input('radio_UEP', 'value'),
    Input('slider_Vh', 'value'),
    Input('radio_embarc', 'value'),
    Input('slider_Hsig', 'value'),
    Input('slider_penTcg', 'value'),
    Input('slider_penTcl', 'value'),
    Input('slider_penEcl', 'value'),
    Input('slider_penMom', 'value'),
    Input('slider_penFht', 'value'),
    Input('slider_penFhd', 'value'),
    Input('slider_penFator', 'value'),
    Input('In_Ecm', 'value'),
    Input('In_El', 'value'),
    Input('In_Ep', 'value'),
    Input('slider_H', 'value'),
    ])

def mostra_modelo(
    combo_modelos,
    combo_ang_raio,
    combo_param,
    radio_mancal,
    radio_UEP,
    slider_Vh,
    radio_embarc,
    slider_Hsig,
    slider_penTcg,
    slider_penTcl,
    slider_penEcl,
    slider_penMom,
    slider_penFht,
    slider_penFhd,
    slider_penFator,
    In_Ecm,
    In_El,
    In_Ep,
    slider_H
    ):
    df = df_dados.set_index('Modelo')
    
    #Associa os valores do banco de dados às variáveis
    teta = array(ast.literal_eval(df.loc[combo_modelos]['teta']))
    Pc = array(ast.literal_eval(df.loc[combo_modelos]['Pc']))
    a = df.loc[combo_modelos]['a']
    b = df.loc[combo_modelos]['b']
    V = df.loc[combo_modelos]['V']
    H = df.loc[combo_modelos]['H']
    G = df.loc[combo_modelos]['G']
    L = df.loc[combo_modelos]['L']
    N = df.loc[combo_modelos]['N']
    S = df.loc[combo_modelos]['S']
    J = df.loc[combo_modelos]['J']
    Pbola = df.loc[combo_modelos]['Pbola']
    Npl = df.loc[combo_modelos]['Npl']
    Nrl = df.loc[combo_modelos]['Nrl']
    Npm = df.loc[combo_modelos]['Npm']
    Nrm = df.loc[combo_modelos]['Nrm']
    Pplat = df.loc[combo_modelos]['Pplat']
    Pcp = df.loc[combo_modelos]['Pcp']
    Dcp = df.loc[combo_modelos]['Dcp']
    Dplat = df.loc[combo_modelos]['Dplat']
    tipocorda = df.loc[combo_modelos]['tipocorda']
    hl = df.loc[combo_modelos]['hl']
    ll = df.loc[combo_modelos]['ll']
    Dc = df.loc[combo_modelos]['Dc']
    t = df.loc[combo_modelos]['t']
    matcordas = df.loc[combo_modelos]['matcordas']
    Tec = df.loc[combo_modelos]['Tec']
    Npend = df.loc[combo_modelos]['Npend']
    Dcabolan = df.loc[combo_modelos]['Dcabolan']
    Dcabopend = df.loc[combo_modelos]['Dcabopend']
    Dcabomoi = df.loc[combo_modelos]['Dcabomoi']
    Pl  = df.loc[combo_modelos]['Pl']
    M  = df.loc[combo_modelos]['M']
    Lpend = df.loc[combo_modelos]['Lpend']
    Ljib = df.loc[combo_modelos]['Ljib']
    Sjib = df.loc[combo_modelos]['Sjib']
    CC1 = df.loc[combo_modelos]['CC1']
    CC2 = df.loc[combo_modelos]['CC2']
    CC3 = df.loc[combo_modelos]['CC3']
    CC4 = df.loc[combo_modelos]['CC4']
    D1 = df.loc[combo_modelos]['D1']
    D2 = df.loc[combo_modelos]['D2']
    D3 = df.loc[combo_modelos]['D3']
    D4 = df.loc[combo_modelos]['D4']
    Elan = df.loc[combo_modelos]['Elan']
    L = df.loc[combo_modelos]['L']
    Pmoi = df.loc[combo_modelos]['Pmoi']
    tetac = radians(df.loc[combo_modelos]['tetac'])

    '''Cálculo do Cvon'''
    #Parâmetros que são calculados de acordo com o Hsig
    Hsig = slider_Hsig
    Vc, Av, CHA, List, Trim  = 0, 0, 0, 0, 0
    if radio_UEP == 'Fixa':
        Vc, Av, CHA, List, Trim = 0, 0, 0, .5, .5
    elif radio_UEP == 'SS':
        Vc = 0.025 * Hsig * Hsig
        Av = 0.0007 * Hsig * Hsig
        if (Av < .07):
            Av = .07
        CHA = 0.007 * Hsig
        if (CHA < .03):
            CHA = .03
        List, Trim = 1.5, 1.5
    elif radio_UEP == 'FPSO':
        Vc = 0.05 * Hsig * Hsig
        Av = 0.0012 * Hsig * Hsig
        if (Av < .07):
            Av = .07
        CHA = .01*(Hsig)**1.1
        if (CHA < .03):
            CHA = .03
        List, Trim = 2.5, 2.5 #Conferir!!!
    
    def func(Pc, teta):

        teta_rad = radians(teta)
        
        #Cálculo do raio
        r = J + L * cos(teta_rad) + S * sin(teta_rad)

        #Tamanho dos vetores
        if type(Pc) is np.ndarray:
            tamanho = len(teta)

        #Raio do jib
        rjib = r + Ljib * cos(teta_rad) + Sjib * sin(teta_rad)
        
        #Parâmetros geométricos
        Dg = (H**2 + (V + G/2)**2)**.5
        tetag = arctan((V + G/2) / H)
        Lcg = (Dg**2 + L**2 - 2 * L * cos(pi - teta_rad - tetag))**.5
        alfa = arcsin((Dg * sin(pi - teta_rad - tetag)) / Lcg)

        #Para guindastes com guincho de lança na lança
        if (H < .01 and V < .01 and G < .01):
            if type(Pc) is np.ndarray:
                for i in range(tamanho):
                    Lcg[i] = L #Considera o guincho de lança no pino do pé
            else:
                Lcg = L
            alfa = 0 * teta_rad #Anula o vetor alfa

        Dl = (a**2 + b**2)**.5
        tetal = arctan(b / a)
        Lcl = (Dl**2 + (L-N)**2 - 2 * Dl * (L-N) * cos(pi - teta_rad - tetal))**.5
        beta = arcsin((Dl * sin(pi - teta_rad - tetal)) / Lcl)

        #Cálculo do Fast Line Factor
        K = 0
        if radio_mancal=='Rolamento':
            K = 1.04
        else:
            K = 1.09
        
        Efm = (K**Npm - 1) / (K**Nrm * Npm * (K - 1))
        Efl = (K**Npl - 1) / (K**Nrl * Npl * (K - 1))

        FLFm = 1 / (Npm * Efm)
        FLFl = 1 / (Npl * Efl)

        alfacl = teta_rad + tetal + tetac - beta - pi/2


        cvon = 1.373 - ((Pc+Pmoi)*2.204623)/(1173913) + Av
        #Limita o valor do Cvon
        if type(Pc) is np.ndarray:
            for i in range(tamanho):
                if (cvon[i] <= (1.1+Av)):
                    cvon[i] = 1.1+Av #adm
                elif (cvon[i] >= (1.33+Av)):
                    cvon[i] = 1.33+Av #adm
        else:
            if (cvon <= (1.1+Av)):
                cvon = 1.1+Av #adm
            elif (cvon >= (1.33+Av)):
                cvon = 1.33+Av #adm

        #Factored Load
        FLkgf = (Pc + Pmoi) * cvon

        #Cálculo do esforço no cabo do moitão
        Tcg = FLkgf * FLFm

        #Cálculo de esforço no cabo da lança
        Esl = (Pl * M * cos(teta_rad) + FLkgf * (L * cos(teta_rad)) + Pbola * (L + Ljib) * cos(teta_rad) + (CC1*D1 + CC2*D2 + CC3*D3) * cos(teta_rad) - Tcg * L * sin(alfa)) / ((L - N) * sin(beta))
        Esl = Esl / Efl
        Tcl = Esl * FLFl

        #Reações no pino do pé da lança
        Rpx = Esl * cos(beta) + Tcg * cos(alfa) + ((CC1+CC2+CC3) + Pl + FLkgf + Pbola) * sin(teta_rad)
        Rpy = (CC1+CC2+CC3+ Pl + FLkgf + Pbola) * cos(teta_rad) - Esl * sin(beta) - Tcg * sin(alfa)
        Rp = (Rpx**2 + Rpy**2)**.5
        gama = arctan(abs(Rpy)/abs(Rpx))

        """Esforço de compressão da lança"""
        Ecl = Rp * cos(gama)

        """Momento"""
        Mom = (J + D1*cos(teta_rad))*CC1 + (J + D2 * cos(teta_rad))*CC2 + (J + D3 * cos(teta_rad)) * CC3 + Pl * (J + M * cos(teta_rad)) + FLkgf * r - Pplat * Dplat - Pcp * Dcp
        
        """Cavalete"""
        Fhd = Esl * cos(alfacl) / sin(tetac)
        Fht = Esl * sin(alfacl) + Fhd * cos(tetac) - Tcl

        result = {}

        result['teta_rad'] = teta_rad
        result['r'] = r
        result['rjib'] = rjib
        result['Dg'] = Dg
        result['tetag'] = tetag
        result['Lcg'] = Lcg
        result['alfa'] = alfa
        result['Dl'] = Dl
        result['tetal'] = tetal
        result['Lcl'] = Lcl
        result['beta'] = beta
        result['K'] = K
        result['Efm'] = Efm
        result['Efl'] = Efl
        result['FLFm'] = FLFm
        result['FLFl'] = FLFl
        result['alfacl'] = alfacl
        result['cvon'] = cvon
        result['FLkgf'] = FLkgf
        result['Tcg'] = Tcg
        result['Esl'] = Esl
        result['Tcl'] = Tcl
        result['Rpx'] = Rpx
        result['Rpy'] = Rpy
        result['Rp'] = Rp
        result['gama'] = gama
        result['Ecl'] = Ecl
        result['Mom'] = Mom
        result['Fhd'] = Fhd
        result['Fht'] = Fht

        return result
    
    teta_rad = func(Pc, teta)['teta_rad']
    r = func(Pc,teta)['r']
    rjib = func(Pc,teta)['rjib']
    Dg = func(Pc,teta)['Dg']
    tetag = func(Pc,teta)['tetag']
    Lcg = func(Pc,teta)['Lcg']
    alfa = func(Pc,teta)['alfa']
    Dl = func(Pc,teta)['Dl']
    tetal = func(Pc,teta)['tetal']
    Lcl = func(Pc,teta)['Lcl']
    beta = func(Pc,teta)['beta']
    K = func(Pc,teta)['K']
    Efm = func(Pc,teta)['Efm']
    Efl = func(Pc,teta)['Efl']
    FLFm = func(Pc,teta)['FLFm']
    alfacl = func(Pc,teta)['alfacl']
    cvon = func(Pc,teta)['cvon']
    FLkgf = func(Pc,teta)['FLkgf']
    Tcg = func(Pc,teta)['Tcg']
    Esl = func(Pc,teta)['Esl']
    Tcl = func(Pc,teta)['Tcl']
    Rpx = func(Pc,teta)['Rpx']
    Rpy = func(Pc,teta)['Rpy']
    Rp = func(Pc,teta)['Rp']
    gama = func(Pc,teta)['gama']
    Ecl = func(Pc,teta)['Ecl']
    Mom = func(Pc,teta)['Mom']
    Fhd = func(Pc,teta)['Fhd']
    Fht = func(Pc,teta)['Fht']

    """Vetores ótimos"""
    Pcot = copy(Pc)
    Tcgot = copy(Tcg)
    Tclot = copy(Tcl)
    Eslot = copy(Esl)
    Eclot = copy(Ecl)
    Momot = copy(Mom)
    Fhdot = copy(Fhd)
    Fhtot = copy(Fht)
    Rpxot = copy(Rpx)
    Rpyot = copy(Rpy)
    Rpot = copy(Rp)
    gamaot = copy(gama)
    FLkgfot = copy(FLkgf)
    cvonot = copy(cvon)

    penTcg = slider_penTcg/100
    penTcl = slider_penTcl/100
    penEcl = slider_penEcl/100
    penMom = slider_penMom/100
    penFht = slider_penFht/100
    penFhd = slider_penFhd/100
    penFator = slider_penFator/100

    """Evita que o programa entre sempre no loop de otimização, o que o deixa lento"""
    if ((slider_penTcg == 100) and (slider_penTcl == 100) and (slider_penEcl == 100) and (slider_penMom == 100) and (slider_penFht == 100) and (slider_penFhd == 100) and (slider_penFator == 100)):
        pass
    
    else:
        """Otimização da tabela"""
        for i in range(len(Pcot)):
            #print(i,(Tcgot[i] > (max(Tcg)*penTcg)) or (Tclot[i] > (max(Tcl)*penTcl)) or (Eclot[i] > (max(Ecl)*penEcl)) or (Momot[i] > (max(Mom)*penMom)) or (Fhtot[i] > (max(Fht)*penFht)) or (Fhdot[i] > (max(Fhd)*penFhd)) or (Pcot[i] > (Pc[i]*penFator)), cont)

            while((Tcgot[i] > (max(Tcg)*penTcg)) or (Tclot[i] > (max(Tcl)*penTcl)) or (Eclot[i] > (max(Ecl)*penEcl)) or (Momot[i] > (max(Mom)*penMom)) or (Fhtot[i] > (max(Fht)*penFht)) or (Fhdot[i] > (max(Fhd)*penFhd)) or (Pcot[i] > (Pc[i]*penFator))):
                Pcot[i] -= 10
                result_ot = func(Pcot[i], teta[i])
                Tcgot[i] = result_ot['Tcg']
                Tclot[i] = result_ot['Tcl']
                Eclot[i] = result_ot['Ecl']
                Momot[i] = result_ot['Mom']
                Fhdot[i] = result_ot['Fhd']
                Fhtot[i] = result_ot['Fht']

    #CÁLCULOS DOS CRITÉRIOS DEFINIDOS NA API 2C
    Vd = 0
    Vc = 0
    #Av = 0
    CHA = 0
    List = 0
    Trim = 0
    
    #Hsig e gravidade em unidades imperiais
    Hsig = Hsig *3.28083 #ft
    grav = 32.2 #ft/s²
    esc_UEP = radio_UEP

    #APLICAÇÃO DAS FUNÇÕES
    esc_embarc = radio_embarc
    Vd = Funcoes.calc_Vd(esc_embarc,Hsig)
    Vc = Funcoes.calc_Vc(esc_UEP,Hsig)
    CHA = Funcoes.calc_CHA(esc_UEP,Hsig)
    List = Funcoes.calc_List(esc_UEP)
    Trim = Funcoes.calc_Trim(esc_UEP)
    Vhmin = Funcoes.calc_Vhmin(Hsig)

    Vh = slider_Vh

    Vr = Vh + (Vd**2+Vc**2)**.5

    #CÁLCULOS DE RIGIDEZ
    El = float(In_El)
    Em = float(In_Ecm)
    Ep = float(In_Ep)
    Alt = slider_H

    Lclan = Lcl - Lpend #m
    kcabolan = (1/9.81) * Npl * El*(1e6) * (.66*.25*(pi*Dcabolan**2)/Lclan) #kgf/m
    kcabopend = (1/9.81) * Npend * Ep*(1e6) * (.66*.25*(pi*Dcabopend**2)/Lpend) #kgf/m
    Htip = Alt +L*sin(teta_rad) #m
    Acb = (.77*.25*(pi*Dcabomoi**2))

    ksustlan = (kcabopend*kcabolan)/(kcabopend+kcabolan) #kgf/m

    Alan = 0

    if (tipocorda == "C"):
        Alan = pi * (Dc**2 - (Dc-2*t)**2) / 4 #m²
        #print(Alan)
    elif (tipocorda == "Q"):
        Alan = Dc**2 - (Dc - t)**2 #m²
    elif (tipocorda == "L"):
        Alan = Dc * 2 * t - t**2 #m²

    eixo_y = {
        '[Gráfico] Ângulo alfa':rad2deg(alfa),
        '[Gráfico] Ângulo beta':rad2deg(beta),
        '[Gráfico] Capacidade Estática':Pc,
        '[Gráfico] Cabo do sistema principal':Tcg,
        '[Gráfico] Cabo de lança':Tcl,
        '[Gráfico] Raio':r,
        '[Gráfico] Sustentação da lança':Esl,
        '[Gráfico] cvon':cvon,
        '[Gráfico] cvonot':cvonot,
        '[Gráfico] Momento':Mom,
        '[Gráfico] Esforço nas hastes traseiras do cav.':Fht,
        '[Gráfico] Esforço nas hastes dianteiras do cav.':Fhd,
        '[Gráfico] Rp':Rp,
        '[Gráfico] Rpx':Rpx,         
        '[Gráfico] Rpy':Rpy,         
        '[Gráfico] alfacl':degrees(alfacl),
        '[Gráfico] Cabo do sistema principal Orig. x Otim.':Tcgot,
        '[Gráfico] Cabo de lança Orig. x Otim.':Tclot,
        '[Gráfico] Sustentação da lança Orig. x Otim.':Tclot,
        '[Gráfico] Momento Orig. x Otim.':Momot,
        '[Gráfico] Esforço nas hastes traseiras do cav. Orig. x Otim.':Fhtot,
        '[Gráfico] Esforço nas hastes dianteiras do cav. Orig. x Otim.':Fhdot,
        '[Gráfico] Esforço de compressão na lança':Ecl,
        '[Gráfico] Esforço de compressão na lança Orig. x Otim.':Eclot,
        '[Gráfico] gama':degrees(gama),
        '[Gráfico] teta':teta,
        '[Gráfico] Capacidade Estática Orig. x Otim.':Pcot,
    }

    resultados = pd.DataFrame({
        'teta':teta,
        'Pc':Pc,
        'Pcot':Pcot
        })
    
    resultados.to_csv('resultado.csv',sep=',')

    eixo_x = teta
    label_x = 'Ângulo da Lança (°)'

    if (combo_ang_raio == 'Raio(m)'):
        eixo_x = r
        label_x = 'Raio (m)'
    
    df.loc[combo_modelos]['[Gráfico] Raio'] = ast.literal_eval(df.loc[combo_modelos]['teta'])

    if combo_param[:9] == '[Gráfico]':

        if combo_param[len(combo_param)-5:] == 'Otim.':

            return dcc.Graph(
                id='grafico_tabela',
                figure={
                    'data':[
                        go.Scatter(
                            name='Original',
                            x=eixo_x,
                            y=eixo_y[combo_param[:len(combo_param)-14]],
                            text='Textoo',
                            mode='markers+lines'
                        ),
                        go.Scatter(
                            name='Otimizado',
                            x=eixo_x,
                            y=eixo_y[combo_param],
                            text='Textoo',
                            mode='markers+lines'
                        ),
                    ],
                    'layout':go.Layout(
                        title=combo_param[10:]+' - '+combo_modelos,
                        xaxis={'title':label_x},
                        yaxis={'title':combo_param[10:]},
                        hovermode='closest'
                    )
                }
            )
        
        else:
            return dcc.Graph(
                id='grafico_tabela',
                figure={
                    'data':[
                        go.Scatter(
                                #name='Otimizado',
                                x=eixo_x,
                                y=eixo_y[combo_param],
                                text='Textoo',
                                mode='markers+lines'
                            ),
                    ],
                    'layout':go.Layout(
                            title=combo_param[10:]+' - '+combo_modelos,
                            xaxis={'title':label_x},
                            yaxis={'title':combo_param[10:]},
                            hovermode='closest'
                    )
                }
            )

    if combo_param[:10] == '[Variável]':

            return str(eval(combo_param[11:]))

@app.callback(
    dash.dependencies.Output('slider_Hsig_out', 'children'),
    [dash.dependencies.Input('slider_Hsig', 'value')])
def update_output(value):
    return 'Hsig = {}m'.format(value)

@app.callback(
    dash.dependencies.Output('slider_Vh_out', 'children'),
    [dash.dependencies.Input('slider_Vh', 'value')])
def update_output2(value):
    return 'Vh = {}ft/s'.format(value)

if __name__ == '__main__':
    app.run_server(debug=True)
