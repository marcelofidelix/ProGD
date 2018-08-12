from numpy import radians, cos, sin, argmax, round_
#import matplotlib.pylab as plt
'''
#Função que gera o gráfico de uma variável
def gera_grafico(linhas,colunas,posicao,titulo,par1,par2):
    plt.subplot(linhas,colunas,posicao)
    plt.plot(par1,par2,lw=3)
    plt.plot(par1,par2,"kx")
    plt.title(titulo,fontsize=10)
    plt.xlim(min(par1)-10,max(par1)*1.1)
    plt.ylim(0,max(par2)*1.1)
    plt.grid()

#Função que gera o gráfico de duas variáveis
def gera_grafico2(linhas,colunas,posicao,titulo,par1,par2,par3):
    plt.subplot(linhas,colunas,posicao)
    plt.plot(par1,par2,lw=3)
    plt.plot(par1,par3,lw=3)
    plt.plot(par1,par3,"kx")
    plt.title(titulo,fontsize=10)
    plt.xlim(min(par1)-10,max(par1)*1.1)
    plt.ylim(0,max(par2)*1.1)
    plt.text(min(par1),max(par2)*.05,'Máximo = '+str(round(max(par3),2)))
    plt.grid()

#Função que gera o gráfico por altura de onda
def graf_Hsig(linhas,colunas,posicao,titulo,par1,par2,par3,par4,par5,par6,par7,par8,par9,par10):
    plt.subplot(linhas,colunas,posicao)
    plt.plot(par1,par2,lw=3,label='Onboard Orig.')
    plt.plot(par1,par3,lw=3,label='Onboard Otm.')
    plt.plot(par1,par3,'kx')
    plt.plot(par1,par10,label='Off Hsig',lw=3)
    plt.plot(par1,par4,label='API 2C 0,5m',lw=.5)
    plt.plot(par1,par5,label='API 2C 1,0m',lw=.5)
    plt.plot(par1,par6,label='API 2C 1,5m',lw=.5)
    plt.plot(par1,par7,label='API 2C 2,0m',lw=.5)
    plt.plot(par1,par8,label='API 2C 2,5m',lw=.5)
    plt.plot(par1,par9,label='API 2C 3,0m',lw=.5)
    plt.title(titulo)
    plt.xlim(min(par1)-10,max(par1)*1.1)
    plt.ylim(0,max(par2)*1.1)
    plt.legend(loc=0,fontsize=7)#, borderaxespad=0.)
    plt.grid()
 '''
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
    if (esc_UEP == 'Botton Supported Structure'):
        Av = 0. #*g
    elif (esc_UEP == 'Ship and Barge in Calm Water'):
        Av = 0.
    elif (esc_UEP == 'TLP'):
        Av = 0.003 * Hsig
        if (Av < .07):
            Av = .07
    elif (esc_UEP == 'Spar'):
        Av = 0.003 * Hsig
        if (Av < .07):
            Av = .07
    elif (esc_UEP == 'SS'):
        Av = 0.0007 * Hsig * Hsig
        if (Av < .07):
            Av = .07
    elif (esc_UEP == 'Drill Ship'):
        Av = 0.0012 * Hsig * Hsig
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
'''
#Gera imagem com gráfico e adiciona ao relatório em html
def gera_graf_html_1(par1,par2,nome,titx):
    plt.plot(par1,par2,lw=3,label=titx)
    plt.plot(par1,par2,'kx')
    plt.ylim(0,max(par2)*1.1)
    plt.title('Rigidez Vertical do Guindaste')
    plt.grid()
    plt.tight_layout()
    plt.legend()
    plt.savefig(nome+'.png')
    plt.gcf().clear()

def gera_graf_html_2(par1,par2,par3,nome,titx,tity):
    plt.plot(par1,par2,lw=3,label=titx)
    plt.plot(par1,par2,'kx')
    plt.plot(par1,par3,lw=3,label=tity)
    plt.plot(par1,par3,'kx')
    plt.ylim(0,max(par2)*1.1)
    plt.title('Máximo Ot.= '+str(int(max(par3))))
    plt.grid()
    plt.tight_layout()
    plt.legend()
    plt.savefig(nome+'.png')
    plt.gcf().clear()

def gera_graf_html_2_2(par1,par2,par3,nome,titx,tity):
    plt.plot(par1,par2,lw=3,label=titx)
    plt.plot(par1,par2,'kx')
    plt.plot(par1,par3,lw=3,label=tity)
    plt.plot(par1,par3,'kx')
    if (max(par2) > max(par3)):
        plt.ylim(0,max(par2)*1.1)
    elif (max(par2) <= max(par3)):
        plt.ylim(0,max(par3)*1.1)
    plt.grid()
    plt.tight_layout()
    plt.legend()
    plt.savefig(nome+'.png')
    plt.gcf().clear()

def gera_graf_html_5(par1,par2,par3,par4,par5,par6,nome,tit2,tit3,tit4,tit5,tit6):
    plt.subplot(3,2,1)
    plt.plot(par1,par2,lw=3,label=tit2)
    plt.plot(par1,par2,'kx')
    plt.grid()
    plt.legend()
    plt.subplot(3,2,2)
    plt.plot(par1,par3,lw=3,label=tit3)
    plt.plot(par1,par3,'kx')
    plt.grid()
    plt.legend()
    plt.subplot(3,2,3)
    plt.plot(par1,par4,lw=3,label=tit4)
    plt.plot(par1,par4,'kx')
    plt.grid()
    plt.legend()
    plt.subplot(3,2,4)
    plt.plot(par1,par5,lw=3,label=tit5)
    plt.plot(par1,par5,'kx')
    plt.grid()
    plt.legend()
    plt.subplot(3,1,3)
    plt.plot(par1,par6,lw=3,label=tit6)
    plt.plot(par1,par6,'kx')
    #plt.ylim(0,max(par2)*1.1)
    plt.grid()
    plt.tight_layout()
    plt.legend()
    plt.savefig(nome+'.png')
    plt.gcf().clear()
'''
'''
def relat(  t1,
            t2,
            modelo,
            FLFl,
            df_dados_GD,
            df_tabela,
            df_geom,
            df_geom_vec,
            df_esfOrig,
            df_rest,
            df_API,
            df_rig,
            df_tab_off,
            df_CS_cabos,
            df_Tabelas_Hsig):

    #Escreve o html
    f_html = open("relatorio.html","w")
    #f_html.write("<object type=\"text/html\" data=\"tab1.html\" width=\"732\" height=\"95\" style=\"overflow:hidden; width: 500px; height: 1050px\"></object>")
    t = '<!DOCTYPE html>'
    t += '<html>'
    t += '<head>'
    t += '<meta charset="UTF-8">'
#    t += '<meta http-equiv="content-type" content="text/html;charset=utf-8">'
#    t += '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">'
    t += '</head>'
    t += '<h1>GDCalc</h1>'
    t += '<h2>Relatório de Análise</h2>'
    t += '<h3>Modelo Selecionado: '+modelo+'</h3>'
    t += '<h4>Tempo de execução: '+str(t2 - t1)+'s</h4>'
    t += '<h3>Parâmetros do Modelo</h3>'
    t += df_dados_GD.to_html(justify='left')
    t += '<h4>Desenho com indicação dos parâmetros geométricos</h4>'
    t += '<img src="Img_01.png" alt="Param" style="width:1240px;height:888px;">'
    t += '<h3>Tabela Original e Otimizada</h3>'
    t += df_tabela.to_html(justify='left')

    t += '<h4>Gráfico das Tabelas Original e Otimizada</h4>'
    gera_graf_html_2(   df_tabela['Âng (°)'].values,
                        df_tabela['Cap (kgf)'].values,
                        df_tabela['Cap Ot. (kgf)'].values,
                        'Img_Tabela','Capacidade estática original (kgf)','Capacidade estática otimizada (kgf)')
    t += '<img src="Img_Tabela.png" alt="Img_Tabela" style="width:640px;height:480px;">'

    t += '<h3>Parâmetros Geométricos</h3>'
    t += df_geom.to_html(justify='left')
    t += '<br>'
    t += df_geom_vec.to_html(justify='left')

    t += '<h4>Ângulos alfa e beta</h4>'
    gera_graf_html_2(   df_tabela['Âng (°)'].values,
                        df_geom_vec['beta (°)'].values,
                        df_geom_vec['alfa (°)'].values,
                        'Img_Alfa_Beta','Alfa (°)','Beta (°)')
    t += '<img src="Img_Alfa_Beta.png" alt="Param" style="width:640px;height:480px;">'

    t += '<h3>Restrições Percentuais</h3>'
    t += df_rest.to_html(justify='left')
    t += '<h3>Tabela de Esforços</h3>'
    t += df_esfOrig.to_html(justify='left')

    t += '<h4>Esforço no Cabo do Guincho</h4>'
    gera_graf_html_2(   df_tabela['Âng (°)'].values,
                        df_esfOrig['Tcg (kgf)'].values,
                        df_esfOrig['Tcg_ot (kgf)'].values,
                        'Img_Tcg','Original (kgf)','Otimizado (kgf)')
    t += '<img src="Img_Tcg.png" alt="Img_Tcg" style="width:640px;height:480px;">'

    t += '<h4>Esforço no Cabo da Lança</h4>'
    gera_graf_html_2(   df_tabela['Âng (°)'].values,
                        df_esfOrig['Tcl (kgf)'].values*FLFl,
                        df_esfOrig['Tcl_ot (kgf)'].values*FLFl,
                        'Img_Tcl','Original (kgf)','Otimizado (kgf)')
    t += '<img src="Img_Tcl.png" alt="Img_Tcl" style="width:640px;height:480px;">'

    t += '<h4>Esforço de Compressão da Lança</h4>'
    gera_graf_html_2(   df_tabela['Âng (°)'].values,
                        df_esfOrig['Ecl (kgf)'].values,
                        df_esfOrig['Ecl_ot (kgf)'].values,
                        'Img_Ecl','Original (kgf)','Otimizado (kgf)')
    t += '<img src="Img_Ecl.png" alt="Img_Ecl" style="width:640px;height:480px;">'

    t += '<h4>Momento Resultante sobre o Pedestal</h4>'
    gera_graf_html_2(   df_tabela['Âng (°)'].values,
                        df_esfOrig['Mom (kgf.m)'].values,
                        df_esfOrig['Mom_ot (kgf.m)'].values,
                        'Img_Mom','Original (kgf.m)','Otimizado (kgf.m)')
    t += '<img src="Img_Mom.png" alt="Img_Mom" style="width:640px;height:480px;">'

    t += '<h4>Esforço nas hastes traseiras do cavalete</h4>'
    gera_graf_html_2(   df_tabela['Âng (°)'].values,
                        df_esfOrig['Fht (kgf)'].values,
                        df_esfOrig['Fht_ot (kgf)'].values,
                        'Img_Fht','Original (kgf)','Otimizado (kgf)')
    t += '<img src="Img_Fht.png" alt="Img_Fht" style="width:640px;height:480px;">'

    t += '<h4>Esforço nas hastes dianteiras do cavalete</h4>'
    gera_graf_html_2(   df_tabela['Âng (°)'].values,
                        df_esfOrig['Fhd (kgf)'].values,
                        df_esfOrig['Fhd_ot (kgf)'].values,
                        'Img_Fhd','Original (kgf)','Otimizado (kgf)')
    t += '<img src="Img_Fhd.png" alt="Img_Fhd" style="width:640px;height:480px;">'

    t += '<h3>Parâmetros API 2C</h3>'
    t += df_API.to_html(justify='left')
    t += '<h3>Cálculo de Rigidez</h3>'
    t += '<br>'
    t += df_rig.to_html(justify='left')

    t += '<h4>Gráfico de Rigidez dos Componentes</h4>'
    gera_graf_html_5(   df_tabela['Âng (°)'].values,
                        df_rig['klan'].values,
                        df_rig['kclan'].values,
                        df_rig['kcabopend'].values,
                        df_rig['ksustlan'].values,
                        df_rig['kcmoi'].values,
                        'Img_Rigidez','Lança','Cabo de Lança','Pendentes','Sust. da Lança','Cabo do Moitão')
    t += '<img src="Img_Rigidez.png" alt="Img_Rigidez" style="width:640px;height:480px;">'

    gera_graf_html_1(   df_tabela['Âng (°)'].values,
                        df_rig['kgd (kgf/mm)'].values,
                        'Img_kgd','Rigidez (kgf/mm)')
    t += '<img src="Img_kgd.png" alt="Img_kgd" style="width:640px;height:480px;">'

    t += '<h3>Coeficientes Dinâmicos e Tabela Offboard (Hsig escolhido)</h3>'
    t += df_tab_off.to_html(justify='left')

    gera_graf_html_2_2(   df_tabela['Âng (°)'].values,
                        df_tab_off['Cvon'].values,
                        df_tab_off['Cvoff'].values,
                        'Img_Cv','Cv Onboard','Cv Offboard')
    t += '<img src="Img_Cv.png" alt="Img_Cv" style="width:640px;height:480px;">'

    gera_graf_html_2_2(   df_tabela['Âng (°)'].values,
                        df_tab_off['Pc (kgf)'].values,
                        df_tab_off['Pcoff (kgf)'].values,
                        'Img_Tab_On_Off','Tabela Onboard (kgf)','Tabela Offboard (kgf)')
    t += '<img src="Img_Tab_On_Off.png" alt="Img_Tab_On_Off" style="width:640px;height:480px;">'

    t += '<h3>Coeficientes de Segurança dos Cabos</h3>'
    t += df_CS_cabos.to_html(justify='left')
    t += '<h3>Tabelas por Altura de Onda</h3>'
    t += df_Tabelas_Hsig.to_html(justify='left')
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += ''
    t += '</body>'
    t += '</html>'
    f_html.write(t)
    f_html.close()
'''