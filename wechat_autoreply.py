#!/usr/bin/env python
#-*- coding: utf-8 -*-

import itchat
to_place=input('Where to go: ')
contact=input('leave contact message or not:')
if set(contact) & set(('t','T','y','Y','是')):
    leave_message=True
else:
    leave_message=False
if to_place is not '':
    send_message='【自动回复】对不起，我去'+to_place+'了，暂时无法回复!'+leave_message*' 如有急事，请联系：15600927092；jerryshen68@hotmail.com'
else:
    send_message='【自动回复】对不起，您联系的用户暂时不在微信旁。如有急事，请联系：15600927092；jerryshen68@hotmail.com'
@itchat.msg_register(itchat.content.TEXT)
def print_content(msg):
    return msg.user.send(send_message)
itchat.auto_login(hotReload=True)#,enableCmdQR = True)
itchat.run()
























