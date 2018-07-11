import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input
import plotly.graph_objs as go
import pandas as pd
import ast
from numpy import radians, cos, array, arctan, arcsin, sin, pi, rad2deg, degrees, copy

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
'[Gráfico] Raio',
'[Gráfico] Cabo do sistema principal',
'[Gráfico] Cabo do sistema principal Orig. x Otim.',
'[Gráfico] Cabo de lança',
'[Gráfico] Cabo de lança Orig. x Otim.',
'[Gráfico] cvon',
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
'[Gráfico] Rpx',
'[Gráfico] Rpy',         
'[Gráfico] gama_teta',
'[Gráfico] Tcgx',         
'[Gráfico] Tcgy',
'[Gráfico] Tclx',
'[Gráfico] Tcly',
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

    html.Div('Fator de limitação'),

    dcc.Slider(
        id='slider_penLim',
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

    html.H4('Altura de onda significativa'),

    dcc.Slider(
        id='slider_Hsig',
        min=0,
        max=5,
        step=0.1,
        value=0,
    ),],style={'width': '80%'}),

    html.Div(id='slider_Hsig_out'),

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

    ],
        style={'width': '100%', 'display': 'inline-block','columnCount': 2}),

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
    Input('slider_Hsig', 'value'),
    Input('slider_penTcg', 'value'),
    Input('slider_penTcl', 'value'),
    Input('slider_penEcl', 'value'),
    Input('slider_penMom', 'value'),
    Input('slider_penFht', 'value'),
    Input('slider_penFhd', 'value'),
    Input('slider_penFator', 'value'),
    Input('slider_penLim', 'value'),
    ])

def mostra_modelo(
    combo_modelos,
    combo_ang_raio,
    combo_param,
    radio_mancal,
    radio_UEP,
    slider_Hsig,
    slider_penTcg,
    slider_penTcl,
    slider_penEcl,
    slider_penMom,
    slider_penFht,
    slider_penFhd,
    slider_penFator,
    slider_penLim):
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
    teta_rad = radians(teta)
    
    #Cálculo do raio
    r = J + L * cos(teta_rad)

    #Tamanho dos vetores
    tamanho = len(teta)

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
    cvon = 1.373 - ((Pc+Pmoi)*2.204623)/(1173913) + Av

    #Limita o valor do Cvon
    for i in range(tamanho):
        if (cvon[i] <= (1.1+Av)):
            cvon[i] = 1.1+Av #adm
        elif (cvon[i] >= (3.33+Av)):
            cvon[i] = 1.33+Av #adm

    #Factored Load
    FLkgf = (Pc + Pmoi) * cvon

    #Parâmetros geométricos
    Dg = ((H)**2+(V+G/2)**2)**.5
    tetag = arctan((V+(G/2))/H)
    Lcg = (Dg**2 + L**2 - 2 * Dg * L * cos (pi - teta_rad - tetag))**.5
    alfa = arcsin((Dg*sin(pi - teta_rad - tetag))/Lcg)

    #Para guindastes com guincho de lança na lança
    if (H < .01 and V < .01 and G < .01):
        for i in range(tamanho):
            Lcg[i] = L #Considera o guincho de lança no pino do pé
        alfa = 0 * teta_rad #Anula o vetor alfa
    Dl = (a**2 + b**2)**.5
    tetal = arctan(b/a)
    Lcl = (Dl**2 + (L-N)**2 - 2 * Dl * (L-N) * cos (pi - teta_rad - tetal))**.5
    beta = arcsin((Dl*sin(pi - teta_rad - tetal))/Lcl)

    #Cálculo do Fast Line Factor
    Kb = 0
    if radio_mancal=='Rolamento':
        Kb = 1.02
    else:
        Kb = 1.045
    FLFl = 1/(Npl*((Kb**Npl - 1) / ((Kb**Nrl) * Npl * (Kb - 1))))
    FLFm = 1/(Npm*((Kb**Npm - 1) / ((Kb**Nrm) * Npm * (Kb - 1))))

    #Cálculo do esforço no cabo do moitão
    Tcg = FLkgf * FLFm

    #Cálculo de esforço no cabo da lança
    Tcl = (Pl*M*cos(teta_rad) + FLkgf*(L*cos(teta_rad) + S*sin(teta_rad)) + Pbola*((L + Ljib)*cos(teta_rad) + Sjib*sin(teta_rad)) + (CC1*D1 + CC2*D2 + CC3*D3 + CC4*D4)*cos(teta_rad) - Tcg*L*sin(alfa)) / ((L - N)*sin(beta))

    #Componentes de Tcg e Tcl 
    Tcgx = -Tcg*(cos(pi/2 - pi/2 - teta_rad - alfa))
    #Tcgy = -Tcg*cos(pi/2 - teta_rad + alfa)
    Tclx = -Tcl*(cos(pi/2 - pi/2 - teta_rad - beta))
    #Tcly = -Tcl*cos(pi/2 - teta_rad + beta)

    Tcgy = Tcg*(sin(pi/2 - pi/2 - teta_rad - alfa))
    Tcly = Tcl*(sin(pi/2 - pi/2 - teta_rad - beta))

    #Reações no pino do pé da lança
    Rpx = Tcgx + Tclx
    Rpy = Tcgy + Tcly - (Pl + FLkgf + Pmoi + Pbola + CC1+CC2+CC3+CC4)
    Rp = (Rpx**2 + Rpy**2)**.5
    gama = arctan(abs(Rpy)/abs(Rpx))
    gama_teta = degrees(gama - radians(teta))

    """Esforço de compressão da lança"""
    Ecl = Rp * cos(gama - teta_rad)

    """Momento"""
    Mom = CC1*(D1*cos(teta_rad)+J) + CC2*(D2*cos(teta_rad)+J) + CC3*(D3*cos(teta_rad)+J) + Pl*(M*cos(teta_rad)+J) + (FLkgf+Pmoi)*(L*cos(teta_rad)+S*sin(teta_rad)+J) + Pbola*((L + Ljib)*cos(teta_rad) + Sjib*sin(teta_rad)) - Pcp*Dcp - Pplat*Dplat
    
    """Cavalete"""
    alfacl = teta_rad - beta
    Fhd = (Tcl*FLFl*Npl*sin(.5*pi+alfacl))/(sin(pi-tetac))
    Fht = (Tcl*FLFl*Npl*sin(.5*pi-tetac+alfacl))/(sin(tetac)) - Tcl*FLFl

    """Vetores ótimos"""
    Pcot = copy(Pc)
    Tcgot = copy(Tcg)
    Tclot = copy(Tcl)
    Tcl_cabo = copy(Tcl * FLFl)
    Tclot_cabo = copy(Tcl_cabo)
    Eclot = copy(Ecl)
    Momot = copy(Mom)
    Fhdot = copy(Fhd)
    Fhtot = copy(Fht)
    Tcgxot = copy(Tcgx)
    Tcgyot = copy(Tcgy)
    Tclxot = copy(Tclx)
    Tclyot = copy(Tcly)
    Rpxot = copy(Rpx)
    Rpyot = copy(Rpy)
    Rpot = copy(Rp)
    gamaot = gama
    FLkgfot = copy(FLkgf)
    cvonot = copy(cvon)
    
    cont = 0

    """Evita que o programa entre sempre no loop de otimização, o que o deixa lento"""
    if ((slider_penTcg == 100) and (slider_penTcl == 100) and (slider_penEcl == 100) and (slider_penMom == 100) and 
    (slider_penFht == 100) and (slider_penFhd == 100) and (slider_penFator == 100) and (slider_penLim == 100)):
        pass
    
    else:

        """Otimização da tabela"""
        for i in range(tamanho):
            
            while((Tcgot[i] > max(Tcg)*slider_penTcg/100) or (Tclot_cabo[i] > max(Tcl_cabo)*slider_penTcl/100)
            or (Eclot[i] > max(Ecl)*slider_penEcl/100) or (Momot[i] > max(Mom)*slider_penMom/100)
            or (Fhtot[i] > max(Fht)*slider_penFht/100) or (Fhdot[i] > max(Fhd)*slider_penFhd/100)
            or (Pcot[i] > Pc[i]*slider_penFator/100) or (Pcot[i] > Pc[i]*slider_penLim/100)):

                Pcot[i] = Pcot[i] - 10
                cvonot[i] = 1.373 - ((Pcot[i]+Pmoi)*2.204623)/(1173913) + Av #adm
                FLkgfot = Pcot[i] * cvonot
                Tcgot[i] = (FLkgfot[i] + Pmoi) * FLFm
                Tclot[i] = (Pl*M*cos(teta_rad[i]) + (FLkgfot[i] + Pmoi)*L*cos(teta_rad[i]) + (CC1*D1 + CC2*D2 + 
                CC3*D3 + CC4*D4)*cos(teta_rad[i]) - Tcgot[i]*L*sin(alfa[i])) / ((L - N)*sin(beta[i]))
                Tclot_cabo[i] = Tclot[i] * FLFl
                Tcgxot[i] = -Tcgot[i]*(cos(teta_rad[i]-alfa[i]))
                Tcgyot[i] = -Tcgot[i]*cos(pi/2-teta_rad[i]+alfa[i])
                Tclxot[i] = -Tclot[i]*(cos(teta_rad[i]-beta[i]))
                Tclyot[i] = -Tclot[i]*cos(pi/2-teta_rad[i]+beta[i])
                Rpxot[i] = - (Tcgxot[i] + Tclxot[i])
                Rpyot[i] = - (Tcgyot[i] + Tclyot[i] - (Pl + FLkgfot[i] + Pmoi + CC1+CC2+CC3+CC4))
                Rpot[i] = (Rpxot[i]**2 + Rpyot[i]**2)**.5
                gamaot[i] = arctan(Rpyot[i]/Rpxot[i])
                Eclot[i] = Rpot[i] * cos(gamaot[i] - teta_rad[i])
                Momot[i] = CC1*(D1*cos(teta_rad[i])+J) + CC2*(D2*cos(teta_rad[i])+J) + CC3*(D3*cos(teta_rad[i])+J) + Pl*(M*cos(teta_rad[i])+J) + (FLkgfot[i]+Pmoi)*(L*cos(teta_rad[i])+ S*sin(teta_rad[i])+J) + Pbola*((L + Ljib)*cos(teta_rad[i])) - Pcp*Dcp - Pplat*Dplat
                Fhdot[i] = (Tclot[i]*FLFl*Npl*sin(.5*pi+alfacl[i]))/(sin(pi-tetac))
                Fhtot[i] = (Tclot[i]*FLFl*Npl*sin(.5*pi-tetac+alfacl[i]))/(sin(tetac)) - Tclot[i]*FLFl
                cont = cont + 1
                print(Pc[i])
                print()


    eixo_y = {
        '[Gráfico] Ângulo alfa':rad2deg(alfa),
        '[Gráfico] Ângulo beta':rad2deg(beta),
        '[Gráfico] Capacidade Estática':Pc,
        '[Gráfico] Cabo do sistema principal':Tcg,
        '[Gráfico] Cabo de lança':Tcl*FLFl,
        '[Gráfico] Raio':r,
        '[Gráfico] Sustentação da lança':Tcl,
        '[Gráfico] cvon':cvon,
        '[Gráfico] Momento':Mom,
        '[Gráfico] Esforço nas hastes traseiras do cav.':Fht,
        '[Gráfico] Esforço nas hastes dianteiras do cav.':Fhd,
        '[Gráfico] Tcgx':Tcgx,
        '[Gráfico] Tcgy':Tcgy,
        '[Gráfico] Tclx':Tclx,
        '[Gráfico] Tcly':Tcly,
        '[Gráfico] Rpx':Rpx,         
        '[Gráfico] Rpy':Rpy,         
        '[Gráfico] gama_teta':gama_teta,
        '[Gráfico] alfacl':degrees(alfacl),
        '[Gráfico] Cabo do sistema principal Orig. x Otim.':Tcgot,
        '[Gráfico] Cabo de lança Orig. x Otim.':Tclot*FLFl,
        '[Gráfico] Sustentação da lança Orig. x Otim.':Tclot,
        '[Gráfico] Momento Orig. x Otim.':Momot,
        '[Gráfico] Esforço nas hastes traseiras do cav. Orig. x Otim.':Fhtot,
        '[Gráfico] Esforço nas hastes dianteiras do cav. Orig. x Otim.':Fhdot,
        '[Gráfico] Esforço de compressão na lança':Ecl,
        '[Gráfico] Esforço de compressão na lança Orig. x Otim.':Eclot,
    }

    eixo_x = teta
    label_x = 'Ângulo da Lança (°)'

    if (combo_ang_raio == 'Raio(m)'):
        eixo_x = r
        label_x = 'Raio (m)'
    
    df.loc[combo_modelos]['[Gráfico] Raio'] = ast.literal_eval(df.loc[combo_modelos]['teta'])
    print(slider_penTcg)
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

if __name__ == '__main__':
    #app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'})
    app.run_server(debug=True)
