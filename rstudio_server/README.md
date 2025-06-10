## 确保系统上的所有软件包都是最新版本
<pre>
sudo dnf update
</pre>

## 确认端口未被占用
<pre>
#确保 RStudio Server 默认端口（8787）未被其他进程占用。使用以下命令检查端口：
sudo netstat -tuln | grep 8787
#如果端口已被占用，可以更改 RStudio Server 使用的端口。在 /etc/rstudio/rserver.conf 中添加以下内容来更改端口：
www-port=8788
</pre>


## 常用rstudio-server命令
<pre>
#启动服务器
systemctl start rstudio-server
#查看服务器状态
systemctl status rstudio-server 
#重启RStudio Server
systemctl restart rstudio-server
#停止服务
systemctl stop rstudio-server
</pre>

## 配置R路径/etc/rstudio/rserver.conf
<pre>rsession-which-r=/staging/software/anaconda3/envs/sesame/bin/R</pre>

## SELinux 可能限制访问，临时关闭测试
<pre>
#关闭 SELinux
setenforce 0
#查看 SELinux 状态
getenforce
</pre> 

## 如果遇到登录或启动 R session 的问题，可以看日志：
<pre>
sudo journalctl -u rstudio-server -f
#或者
cat ~/.local/share/rstudio/log/rsession-ruser.log | tail -n 20
</pre>

## 设置端口防火墙
<pre>
#检查防火墙状态
systemctl status firewalld
#若未运行，可先启动
sudo systemctl start firewalld
sudo systemctl enable firewalld
#开放 8787 端口
firewall-cmd --permanent --add-port=8787/tcp
firewall-cmd --reload
#确认端口是否开放
firewall-cmd --list-ports
</pre>

## SSSD (System Security Services Daemon)，它是一个用于管理与远程身份验证服务（如 LDAP 或 Kerberos）交互的守护进程
<pre>
#清除缓存文件
rm -rf /var/lib/sss/db/*
#然后重启 SSSD 服务
systemctl restart sssd
#检查服务状态
systemctl status sssd
#升级 SSSD
dnf upgrade sssd
</pre>