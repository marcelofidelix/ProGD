#coding: utf-8
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input
import plotly.graph_objs as go
import pandas as pd
import ast
from numpy import radians, cos, array, arctan, arcsin, sin, pi, rad2deg, degrees, copy, clip, arccos, round
import numpy as np
import Funcoes
import pdb
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
'[Gráfico] Capacidade onboard x offboard',
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
'[Gráfico] Momento de inércia',
'[Gráfico] Momento de inércia Orig. x Otim.',
'[Gráfico] Sustentação da lança',
'[Gráfico] Sustentação da lança Orig. x Otim.',
'[Gráfico] alfacl',
'[Gráfico] Rp',
'[Gráfico] Rpx',
'[Gráfico] Rpy',
'[Gráfico] gama',
'[Gráfico] teta',
'[Gráfico] Lcg',
'[Gráfico] Lcl',
'[Gráfico] Lclan',
'[Gráfico] kcabolan',
'[Gráfico] kcabopend',
'[Gráfico] Htip',
'[Gráfico] ksustlan',
'[Gráfico] kcabomoi',
'[Gráfico] dsustlan',
'[Gráfico] dlan',
'[Gráfico] dmoi',
'[Gráfico] tetad',
'[Gráfico] dv_sl',
'[Gráfico] dv_lan',
'[Gráfico] dv_moi',
'[Gráfico] dv',
'[Gráfico] kgd',
'[Gráfico] cvoff',
'[Gráfico] eps',
'[Gráfico] Deps',
'[Gráfico] Factored Load',
'[Gráfico] Design Factors (Running Rigging)',
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
'[Variável] slider_penTcg',
'[Variável] slider_penTcl',
'[Variável] slider_penEcl',
'[Variável] slider_penMom',
'[Variável] slider_penFht',
'[Variável] slider_penFhd',
'[Variável] slider_penI',]
#Ordena a lista
lista_param = sorted(lista_param)
#Função que gera tabela a partir de um DataFrame
def generate_table(dataframe, max_rows=100):
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +
        # Body
        [html.Tr([
            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
        ]) for i in range(min(len(dataframe), max_rows))])
#lista de dicionários que vai receber os valores que populam o dropdown
combo_param = []
for i in lista_param:
    combo_param.append({'label':i,'value':i})
#Lista de opções que irão habitar o radio com opções para o eixo x
combo_grafico = [
    {'label': 'Ângulo(°)', 'value': 'Ângulo(°)'},
    {'label': 'Raio(m)', 'value': 'Raio(m)'},
    {'label': 'Carga(kgf)', 'value':'Carga(kgf)'}
]
#Definição do layout do app
app.layout = html.Div([
    #html.Div(children=[html.H4(children='US Agriculture Exports (2011)'),generate_table(df_dados)]),
    html.Div([
    html.H3('GDCalc'),

    html.H4('Equipamento'),

    dcc.Dropdown(
    id='combo_modelos',
    options=combo_modelos_valores,
    value='[UN-BC] American 5750R 90ft'
    ),

    ],style={'width': '45%', 'display': 'inline-block'}),

    html.Div([

    html.Div(
    
    [html.H4('Penalizações'),

    html.Div('Cabo Moitão'),

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

    html.Div('Cabo Lança'),

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

    html.Div('Comp. Lança'),

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

    html.Div('H. tras. cav.'),

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

    html.Div('H. dian. cav.'),

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

    html.Div('M. de I. de giro'),

    dcc.Slider(
        id='slider_penI',
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

    html.Div('F. de pen.'),

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
        value='SS',
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
    Input('slider_penI', 'value'),
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
    slider_penI,
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
    ########################
    #FUNÇÃO COM OS CÁLCULOS#
    ########################
    def func(Pc, teta):
        teta_rad = radians(teta)
        
        #Cálculo do raio
        r = J + L * cos(teta_rad) + S * sin(teta_rad)

        #Cálculo do raio do CG da lança
        rcg = J + M * cos(teta_rad)

        #Tamanho dos vetores
        if type(Pc) is np.ndarray:
            tamanho = len(teta)

        #Raio do jib
        rjib = r + Ljib * cos(teta_rad) + Sjib * sin(teta_rad)
        
        #Parâmetros geométricos
        Dg = (H**2 + (V + G/2)**2)**.5
        try:
            tetag = arctan((V + G/2) / H)
        except:
            tetag = arctan((V + G/2) / 0.1)
        Lcg = (Dg**2 + L**2 - 2 * L * Dg * cos(pi - teta_rad - tetag))**.5
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
            #K = 1.02
        else:
            K = 1.09
            #K = 1.045
        
        Efm = (K**Npm - 1) / (K**Nrm * Npm * (K - 1))
        Efl = (K**Npl - 1) / (K**Nrl * Npl * (K - 1))
        FLFm = 1 / (Npm * Efm)
        FLFl = 1 / (Npl * Efl)

        #Apesar do tetac ser um dado de entrada, vamos redefinir o seu valor!
        tetac = arctan(a/b)

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
        #Esl = Esl / Efl
        Tcl = Esl / (Npl * Efl)
        #Uma vez calculado Tcl, podemos aplicar o o FLFl médio no Esl
        Esl = Esl * (((1/Efl - 1)/2) + 1)
        #Reações no pino do pé da lança
        Rpx = Esl * cos(beta) + Tcg * cos(alfa) + ((CC1+CC2+CC3) + Pl + FLkgf + Pbola) * sin(teta_rad)
        Rpy = (CC1+CC2+CC3+ Pl + FLkgf + Pbola) * cos(teta_rad) - Esl * sin(beta) - Tcg * sin(alfa)
        #Gambiarra para considerar o Clarke Chapman, que possui roldanas para o moitão no cavalete!
        if max(alfa) > 15.5:
            Rpx = Esl * cos(beta) + 4*FLFm*Tcg * cos(alfa) + ((CC1+CC2+CC3) + Pl + FLkgf + Pbola) * sin(teta_rad)
            Rpy = (CC1+CC2+CC3+ Pl + FLkgf + Pbola) * cos(teta_rad) - Esl * sin(beta) - 4*FLFm*Tcg * sin(alfa)  
        Rp = (Rpx**2 + Rpy**2)**.5
        gama = abs(arctan(abs(Rpy)/abs(Rpx)) - teta_rad)
        #Cálculo do esforço de compressão aplicado sobre a lança
        Ecl = Rp * cos(gama)
        #Cálculo do momento resultante sobre o pedestal
        Mom = (J + D1*cos(teta_rad))*CC1 + (J + D2 * cos(teta_rad))*CC2 + (J + D3 * cos(teta_rad)) * CC3 + Pl * (J + M * cos(teta_rad)) + FLkgf * r + Pbola * rjib - Pplat * Dplat - Pcp * Dcp
        #Cálculo dos esforços nas hastes do cavalete
        Fhd = Esl * cos(alfacl) / sin(tetac)
        Fht = Esl * sin(alfacl) + Fhd * cos(tetac) - Tcl
        #Cálculo do momento de inércia de giro do guindaste
        I = Pcp*Dcp**2 + Pplat*Dplat**2 + CC1*(J+D1*cos(teta_rad))**2 + CC2*(J+D2*cos(teta_rad))**2 + CC3*(J+D3*cos(teta_rad))**2 + Pl*(J+rcg)**2 + FLkgf*r**2 + Pbola*rjib**2
        #Dicionário com os resultados de interesse
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
        result['I'] = I

        return result
    
    result = func(Pc, teta)
    
    teta_rad = result['teta_rad']
    r = result['r']
    rjib = result['rjib']
    Dg = result['Dg']
    tetag = result['tetag']
    Lcg = result['Lcg']
    alfa = result['alfa']
    Dl = result['Dl']
    tetal = result['tetal']
    Lcl = result['Lcl']
    beta = result['beta']
    K = result['K']
    Efm = result['Efm']
    Efl = result['Efl']
    FLFm = result['FLFm']
    FLFl = result['FLFl']
    alfacl = result['alfacl']
    cvon = result['cvon']
    FLkgf = result['FLkgf']
    Tcg = result['Tcg']
    Esl = result['Esl']
    Tcl = result['Tcl']
    Rpx = result['Rpx']
    Rpy = result['Rpy']
    Rp = result['Rp']
    gama = result['gama']
    Ecl = result['Ecl']
    Mom = result['Mom']
    Fhd = result['Fhd']
    Fht = result['Fht']
    I = result['I']

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
    Iot = copy(I)
    gamaot = copy(gama)
    FLkgfot = copy(FLkgf)
    cvonot = copy(cvon)

    penTcg = slider_penTcg/100
    penTcl = slider_penTcl/100
    penEcl = slider_penEcl/100
    penMom = slider_penMom/100
    penFht = slider_penFht/100
    penFhd = slider_penFhd/100
    penI = slider_penI/100
    penFator = slider_penFator/100
    """Evita que o programa entre sempre no loop de otimização, o que o deixa lento"""
    if ((slider_penTcg == 100) and (slider_penTcl == 100) and (slider_penEcl == 100) and (slider_penMom == 100) and 
        (slider_penFht == 100) and (slider_penFhd == 100) and (slider_penFator == 100) and (slider_penI == 100)):
        pass
    else:
        """Otimização da tabela"""
        for i in range(len(Pcot)):
            while((Tcgot[i] > (max(Tcg)*penTcg)) or (Tclot[i] > (max(Tcl)*penTcl)) or (Eclot[i] > (max(Ecl)*penEcl)) or 
                (Momot[i] > (max(Mom)*penMom)) or (Fhtot[i] > (max(Fht)*penFht)) or (Fhdot[i] > (max(Fhd)*penFhd)) or 
                (Pcot[i] > (Pc[i]*penFator)) or (Iot[i] > max(I)*penI)):
                Pcot[i] -= 10
                result_ot = func(Pcot[i], teta[i])
                Tcgot[i] = result_ot['Tcg']
                Tclot[i] = result_ot['Tcl']
                Eclot[i] = result_ot['Ecl']
                Momot[i] = result_ot['Mom']
                Fhdot[i] = result_ot['Fhd']
                Fhtot[i] = result_ot['Fht']
                Iot[i] = result_ot['I']
    ############################################
    #CÁLCULOS DOS CRITÉRIOS DEFINIDOS NA API 2C#
    ############################################
    Vd = 0
    Vc = 0
    #Av = 0
    CHA = 0
    List = 0
    Trim = 0
    #Hsig e gravidade em unidades imperiais#
    Hsig = Hsig *3.28083 #ft
    grav = 32.2 #ft/s²
    esc_UEP = radio_UEP
    #######################
    #APLICAÇÃO DAS FUNÇÕES#
    #######################
    esc_embarc = radio_embarc
    Vd = Funcoes.calc_Vd(esc_embarc,Hsig)
    Vc = Funcoes.calc_Vc(esc_UEP,Hsig)
    CHA = Funcoes.calc_CHA(esc_UEP,Hsig)
    List = Funcoes.calc_List(esc_UEP)
    Trim = Funcoes.calc_Trim(esc_UEP)
    Vhmin = Funcoes.calc_Vhmin(Hsig)
    Vh = slider_Vh
    Vr = Vh + (Vd**2+Vc**2)**.5
    #####################
    #CÁLCULOS DE RIGIDEZ#
    #####################
    #Módulos de elasticidade dos cabos, conforme entrada na GUI
    El = float(In_El) * 10**6 #Pa
    Em = float(In_Ecm) * 10**6 #Pa
    Ep = float(In_Ep) * 10**6 #Pa
    #Distância vertical entre o pino do pé da lança e o mar
    Alt = slider_H
    #Comprimento entre as selas fixa e flutuante medido ao longo dos cabos de lança
    Lclan = Lcl - Lpend #m
    #Rigidez dos cabos de lança
    #Considera o fator de enchimento de 0,68 do cabo Briden Dyform 8 8x26
    kcabolan = (1/9.81) * Npl * El * (.68*.25*(pi*Dcabolan**2)/Lclan) #kgf/m
    #Rigidez dos pendentes
    kcabopend = (1/9.81) * Npend * Ep * (.68*.25*(pi*Dcabopend**2)/Lpend) #kgf/m
    #Distância vertical entre a ponta da lança e a água
    Htip = Alt + L*sin(teta_rad) - S*sin(teta_rad)#m
    #Área da seção do cabo do moitão
    #Considera o fator de enchimento de 0.77 do cabo Bridon Dyform 34LR 34x7
    Acb = (.77*.25*(pi*Dcabomoi**2))
    #Rigidez do sistema de sustentação da lança
    ksustlan = (kcabopend*kcabolan)/(kcabopend+kcabolan) #kgf/m
    #Cálculo da área de uma corda
    Alan = 0
    if (tipocorda == "C"):
        Alan = pi * (Dc**2 - (Dc-2*t)**2) / 4 #m²
    elif (tipocorda == "Q"):
        Alan = Dc**2 - (Dc - t)**2 #m²
    elif (tipocorda == "L"):
        Alan = Dc * 2 * t - t**2 #m²
    #Rigidez das 4 cordas
    klan = (1/9.81) * ((200e9) * 4 * Alan) / (L) #kgf/m
    #Cálculo da rigidez dos cabos do moitão no trecho entre a ponta de a lança e o moitão
    kcabomoi_ae = (1/9.81) * Npm * Em * (Acb / Htip) #kgf/m
    #Cálculo da rigidez dos cabos do moitão no trecho entre o guincho e a ponta da lança
    kcabomoi_guin = (1/9.81) * Em * (Acb / Lcg) #kgf/m
    #Cálculo da rigidez total dos cabos do moitão (associação em série)
    kcabomoi = (kcabomoi_ae*kcabomoi_guin)/(kcabomoi_ae+kcabomoi_guin)
    ###########################
    #CÁLCULO DOS DESLOCAMENTOS#
    ###########################
    #Deslocamento total do sistema de sustentação da lança (Cabo de lança e pendentes)
    dsustlan = Tclot / ksustlan #m
    #Deslocamento da lança sob compressão
    dlan = Eclot / klan #m
    #Deslocamento total do cabo do moitão
    dmoi = (Pcot + Pmoi) / kcabomoi #m
    #Deslocamento angular da lança devido a todas as deformações combinadas
    #tetad = pi - tetal - arccos((Dl**2 + (L-N)**2 - (Lcl * dsustlan)**2) / (2*Dl*(L-N))) #rad
    eps = pi - tetal - teta_rad
    #Variação do ângulo da lança devido à extensão do sistema de sust. da lança
    Deps = arccos(((Lcl+dsustlan)**2 - Dl**2 - (L - N)**2) / (-2*Dl*(L - N))) - eps
    #Ângulo da lança após a extensão do sistema de sustentação da lança
    tetad = teta_rad - Deps
    #Deslocamento vertical da ponta da lança
    #dv_sl = L*(sin(teta_rad) - sin(teta_rad - Deps))
    dv_sl = L * (teta_rad - tetad) * cos(teta_rad) #Devido à extensão do sistema de lança
    dv_lan = dlan * sin(teta_rad) #Devido à compressão da lança
    dv_moi = dmoi #Devido à extensão dos cabos do moitão
    dv = dv_sl + dv_lan + dv_moi
    #Rigidez do guindaste
    kgd = (Pcot + Pmoi) / dv #kgf/m
    #Rigidez em lbf/ft
    kgd = kgd * 0.671968975 #lbf/ft
    ########################
    #CÁLCULO DO Cv OFFBOARD#
    ########################
    #Cálculo do Cv offboard
    cvoff = 1 + Vr * ((kgd)/(grav*((Pcot+Pmoi)*2.204623)))**.5 #adm
    #Relação entre os coeficientes dinâmicos on e offboard
    rel = cvonot/cvoff #amd
    for i in range(len(Pcot)):
        if (rel[i] > 1):
            rel[i] = 1
    #Cálculo da tabela de carga offboard
    Pcoff = Pcot * rel #kgf

    #Coeficiente de Segurança dos Cabos
    CS_cabos = 10000 / (0.004 * (Pc * 2.20462262) + 1910)
    
    for i in range(len(Pcot)):
        if CS_cabos[i] < 3:
            CS_cabos[i] = 3
        if CS_cabos[i] > 5:
            CS_cabos[i] = 5

    eixo_y = {
        '[Gráfico] Ângulo alfa':rad2deg(alfa),
        '[Gráfico] Ângulo beta':rad2deg(beta),
        '[Gráfico] Capacidade Estática':Pc,
        '[Gráfico] Capacidade Estática Orig. x Otim.':Pcot,
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
        '[Gráfico] tetad':tetad,
        '[Gráfico] Momento de inércia':I,
        '[Gráfico] Momento de inércia Orig. x Otim.':Iot,
        '[Gráfico] Lcl':Lcl,
        '[Gráfico] Lcg':Lcg,
        '[Gráfico] Lclan':Lclan,
        '[Gráfico] kcabolan':kcabolan,
        '[Gráfico] Htip':Htip,
        '[Gráfico] ksustlan':ksustlan,
        '[Gráfico] kcabomoi':kcabomoi,
        '[Gráfico] dsustlan':dsustlan,
        '[Gráfico] dlan':dlan,
        '[Gráfico] dmoi':dmoi,
        '[Gráfico] dv_sl':dv_sl,
        '[Gráfico] dv_lan':dv_lan,
        '[Gráfico] dv_moi':dv_moi,
        '[Gráfico] dv':dv,
        '[Gráfico] kgd':kgd,
        '[Gráfico] eps':degrees(eps),
        '[Gráfico] Deps':degrees(Deps),
        '[Gráfico] cvoff':cvoff,
        '[Gráfico] Capacidade onboard x offboard':Pcoff,
        '[Gráfico] Factored Load':Pcot*cvon,
        '[Gráfico] Design Factors (Running Rigging)': CS_cabos
    }
    resultados = pd.DataFrame({
        'teta':teta,
        'raio':r,
        'Pc':round(Pc),
        'Pcot':round(Pcot),
        'Pcoff':round(Pcoff)
        })

    resultados.to_csv('resultado.csv', sep=',')

    eixo_x = teta
    label_x = 'Ângulo da Lança (°)'

    if (combo_ang_raio == 'Raio(m)'):
        eixo_x = r
        label_x = 'Raio (m)'
    if (combo_ang_raio == 'Carga(kgf)'):
        eixo_x = Pcot
        label_x = 'Carga no gancho (kgf)'
    
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
                        plot_bgcolor='rgb(229,229,229)',
                        title=combo_param[10:]+' - '+combo_modelos,
                        xaxis={'title':label_x},
                        yaxis={'title':combo_param[10:]},
                        hovermode='closest'
                    )
                }
            )
        
        if combo_param[len(combo_param)-4:] == 'oard':

            return dcc.Graph(
                id='grafico_tabela',
                figure={
                    'data':[
                        go.Scatter(
                            name='Onboard',
                            x=eixo_x,
                            y=eixo_y['[Gráfico] Capacidade Estática Orig. x Otim.'],
                            text='Textoo',
                            mode='markers+lines'
                        ),
                        go.Scatter(
                            name='Offboard',
                            x=eixo_x,
                            y=eixo_y['[Gráfico] Capacidade onboard x offboard'],
                            text='Textoo',
                            mode='markers+lines'
                        ),
                    ],
                    'layout':go.Layout(
                        plot_bgcolor='rgb(229,229,229)',
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
                            plot_bgcolor='rgb(229,229,229)',
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
