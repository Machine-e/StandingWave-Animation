import matplotlib.pyplot as plt
from matplotlib.widgets import Button,Slider
import numpy as np
import matplotlib.animation as animation

#不考虑泛音，只计算基频

#参数初始化
now_wall_x = 100     #墙的位置
now_source_x = 0   #源点的位置
now_t = 0    #时间

A = 5  # 振幅
T = 1000      #张力

air_density = 1.225    #空气密度 (kg/m^3)
gama = 1.4 #干燥空气的比热容比(定压比热容/定容比热容)
R =  287 #特定气体常数(J/(kg·K))
Temp = 293 #20℃的开尔文温度(K)

u = ((gama*R*Temp) / air_density) ** (1/2) #波速 (m/s)

string_real_density_steel = 8900  #钢弦的密度 (kg/m^3)
string_real_density_nylon = 1300  #尼龙弦的密度 (kg/m^3)

string_density_steel = 0.00586  #钢弦的线密度 (kg/m)
string_density_nylon = 0.00400  #尼龙弦的线密度 (kg/m)
string_density_judge = True  # 弦的种类，1钢，0尼龙，初始为钢弦

    
string_width = 0.001143  #弦的粗细(m)  用的是5弦的直径 音高A4
string_length = 0.648  #弦的长度(m)
string_width_r = string_width / 2  #弦的半径(m)

string_linear_density_steel = string_real_density_steel * ((string_width_r**2 )*np.pi)     #钢弦的线密度 (kg/m)
string_linear_density_nylon = string_real_density_nylon * ((string_width_r**2 )*np.pi)    #尼龙弦的线密度 (kg/m)

if string_density_judge:
    string_density = string_linear_density_steel
else:
    string_density = string_linear_density_nylon
#波疏介质与波密介质的选择
wall_type = True  # 墙的种类，0没有半波损失，1有，初始为波密介质

now_start = False
f = u/(2*string_length)
lambda1 = u/f
w = 2*np.pi*f

class Wave:
    #属性分配
    def __init__(self, source_x, wall_x):
        self.source_x = source_x
        self.wall_x = wall_x
        self.x = np.linspace(source_x, source_x, 1000)
        self.y = None
        self.x_fan = np.linspace(wall_x, wall_x, 1000)
        self.y_fan = None
    #计算y值
    def cal_y(self, t):
        now_x = self.source_x + u * t
        t2 = 0
        have_fan = False
        #波反弹的判断
        if now_x >= self.wall_x:
            now_x = self.wall_x
            t2 = t - (self.wall_x - self.source_x) / u
            have_fan = True
        x_fan = 0
        #x的计算
        if have_fan:
            x_fan = self.wall_x - t2 * u
            if x_fan <= self.source_x:
                x_fan = self.source_x
                self.x_fan = np.linspace(self.source_x, self.wall_x, 1000)
            else:
                self.x_fan = np.linspace(self.wall_x - t2 * u, self.wall_x, 1000)
        self.x = np.linspace(self.source_x, now_x, 1000)
        #反弹前y的计算
        if not have_fan:
            self.y = A * np.cos(w * (t - self.x / u) )
        #反弹后y的计算
        else:
            wave_y = A * np.cos(w * (t - self.x / u) )
            wave_y_fan = None
        #考虑半波损失(dense material)
            if wall_type:
                wave_y_fan = np.where(self.x >= x_fan,
                                      A * np.cos(w * (t - (2 * self.wall_x - self.x) / u) ), 0)
        #不考虑半波损失(rare material)
            else:
                wave_y_fan = np.where(self.x >= x_fan,
                                      A * np.cos(w * (t - (2 * self.wall_x - self.x) / u) + np.pi), 0)
            self.y = wave_y_fan + wave_y


fig, ax = plt.subplots(figsize=(15, 12)) #初始化窗口大小

#定义滑块与按钮
sfreq = None
stension = None
sA = None
sstring_width = None
button1 = None
sair_density = None
su =None #该参数改变波速，但仅用于演示波的轨迹(因为初始波速过快，无法展现出入射波与反射波形成驻波的过程)
def init():
    global sfreq
    global stension
    global sA
    global sstring_width
    global sair_density
    global su

    # 滑块设置
    axfreq = plt.axes([0.25, 0.20, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    sfreq = Slider(axfreq, 'Frequency', 0, 1000, valinit=f)
    sfreq.on_changed(on_sfreq_change)
    
    axtension = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    stension = Slider(axtension, 'Tension', 700,2000, valinit=T)
    stension.on_changed(on_stension_change)

    axA = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    sA = Slider(axA, 'A', 0, 20, valinit=A)
    sA.on_changed(on_sA_change)
    
    axstring_width = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    sstring_width = Slider(axstring_width, 'string_width', 0, 0.01, valinit=string_width)
    sstring_width.on_changed(on_sstring_width_change)
   
    axair_density = plt.axes([0.25, 0.00, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    sair_density = Slider(axair_density, 'air_density', 0.5,50, valinit=air_density)
    sair_density.on_changed(on_sair_density_change)
    #该参数改变波速，但仅用于演示波的轨迹(初始波速过快，无法展现出入射波与反射波形成驻波的过程)
    axu = plt.axes([0.45, 0.95, 0.5, 0.03], facecolor='lightgoldenrodyellow')
    su = Slider(axu, 'wave speed(not object of research,only using to demonstrate the wave trajectory)', 0,400, valinit=u)
    su.on_changed(on_su_change)
    #按钮设置
    button1_ax = plt.axes((0.05,0.20,0.07,0.05))  # 按钮的位置和大小
    button1 = Button(button1_ax, 'Restart')
    button1.on_clicked(on_button1_clicked)
    
    button2_ax = plt.axes((0.10,0.15, 0.02, 0.02))  # 按钮的位置和大小
    button2 = Button(button2_ax, '●')
    button2.on_clicked(on_button2_clicked)
    
    button3_ax = plt.axes((0.10,0.10, 0.02, 0.02)) # 按钮的位置和大小
    button3 =Button(button3_ax, '●')
    button3.on_clicked(on_button3_clicked)
    
    draw_wave(None)
    plt.show()


wave = Wave(now_source_x, now_wall_x)
num = 0

#绘制图像
def draw_wave(frame):
    global now_t
    global num
    plt.subplots_adjust(bottom=0.30)  # 为滑块留出空间
    #设置坐标轴
    ax.clear()
    ax.set_xlim(0,100)
    ax.set_ylim(-15,15)
    ax.set_xlabel('x')
    ax.set_ylabel("y")
    plt.xticks(range(0,100,1))
    plt.yticks(range(-15,15,1))
    #绘图
    wave.cal_y(now_t)
    ax.plot(wave.x, wave.y)
    #标注起止线
    if wall_type:
        ax.plot([now_wall_x, now_wall_x], [-15, 15], color='red', linewidth=2)
    else:
        ax.plot([now_wall_x, now_wall_x], [-15, 15], color='black', linewidth=2)
    ax.plot([now_source_x, now_source_x], [-15, 15], color='green', linewidth=2)
    # 添加带有框的文本框
    if wall_type:
        ax.annotate(" " * (9 - (len("dense") - 1) * 2) + "dense", (0, 0), xytext=(-10,-22.5),
                    bbox=dict(boxstyle='square,pad=0.3', fc='lightgrey', alpha=1, linewidth=0.7), fontsize=10)
    else:
        ax.annotate(" " * (9 - (len("rare") - 1) * 2) + "rare", (0, 0), xytext=(-10,-22.5),
                    bbox=dict(boxstyle='square,pad=0.3', fc='lightgrey', alpha=1, linewidth=0.7), fontsize=10)
    if string_density_judge:
        ax.annotate(" " * (9 - (len("steel") - 1) * 2) + "steel", (0, 0), xytext=(-10,-24.8),
                    bbox=dict(boxstyle='square,pad=0.3', fc='lightgrey', alpha=1, linewidth=0.7), fontsize=10)
    else:
        ax.annotate(" " * (9 - (len("nylon") - 1) * 2) + "nylon", (0, 0), xytext=(-10,-24.8),
                    bbox=dict(boxstyle='square,pad=0.3', fc='lightgrey', alpha=1, linewidth=0.7), fontsize=10)
     #波长显示
    ax.text(0.09,-0.10, 'lambda: {:.5f} m'.format(lambda1), transform=ax.transAxes, fontsize=12,
          )

    num += 1
    now_t += 0.04

def on_sfreq_change(val):
    global now_t
    global f
    global w
    global lambda1
    now_t = 0
    f = sfreq.val
    lambda1 = u/f
    w = 2 * np.pi * f
    
def on_stension_change(val):
    global now_t
    global f
    global w
    global T
    global string_length
    global string_density
    global string_density_judge
    global lambda1
    if string_density_judge:
        string_density = string_linear_density_steel
    else:
        string_density = string_linear_density_nylon
    now_t = 0
    T = stension.val
    f = 1/(2*string_length) *((T/string_density_steel)**(1/2))
    w = 2 * np.pi * f
    lambda1 = u/f
    sfreq.set_val(f)
    
def on_sA_change(val):
    global A
    global now_t
    global f
    global w
    A = sA.val
    now_t = 0
    
def on_sstring_width_change(val):
    global string_width
    global now_t
    global f
    global w
    global string_density
    global lambda1
    string_width = sstring_width.val
    string_width_r = string_width / 2 
    string_linear_density_steel = string_real_density_steel * ((string_width_r**2 )*np.pi)     #钢弦的线密度 (kg/m)
    string_linear_density_nylon = string_real_density_nylon * ((string_width_r**2 )*np.pi)    #尼龙弦的线密度 (kg/m)
    
    if string_density_judge:
        string_density = string_linear_density_steel
    else:
        string_density = string_linear_density_nylon
    
    now_t = 0
    f = 1/(2*string_length) *((T/string_density)**(1/2))
    w = 2 * np.pi * f
    lambda1 = u/f
    sfreq.set_val(f)
    
def on_sair_density_change(val):
    global air_density
    global now_t
    global f
    global w
    global string_density
    global lambda1
    global u
    air_density = sair_density.val
    now_t = 0
    u = ((gama*R*Temp) / air_density) ** (1/2)
    f = u/(2*string_length) 
    lambda1 = u/f
    w = 2 * np.pi * f
    sfreq.set_val(f)
    su.set_val(u)
    
def on_su_change(val):
    global now_t
    global f
    global w
    global string_density
    global lambda1
    global u
    u = su.val
    now_t = 0
    f = u/(2*string_length)
    lambda1 = u/f
    w = 2 * np.pi * f
    sfreq.set_val(f)
    
def on_button1_clicked(event):
    global now_t
    now_t = 0
    global now_start
    now_start = True
def on_button2_clicked(event):
    global now_t
    global wall_type
    now_t = 0
    wall_type = not wall_type
def on_button3_clicked(event):
    global now_t
    global string_density
    global string_density_steel
    global string_density_nylon
    global string_density_judge
    global lambda1

    now_t = 0
    string_density_judge = not string_density_judge
    if string_density_judge:
        string_density = string_linear_density_steel
    else:
        string_density = string_linear_density_nylon

    f = 1/(2*string_length) *((T/string_density)**(1/2))
    w = 2 * np.pi * f
    lambda1 = u/f
    sfreq.set_val(f)
def main():
    global ani
    ani = animation.FuncAnimation(fig, draw_wave, frames=100, interval=20,blit = False)
    init()
    

if __name__ == '__main__':
    main()
