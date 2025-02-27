# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Starterator Dialog Windows

import MySQLdb
from gi.repository import Gtk, GObject, Gdk
import configparser
import os
import subprocess
import sys
import threading
import time
from . import utils

class StarteratorExceptionDialog(Gtk.Dialog):
    def __init__(self, parent):
        Gtk.Dialog.__init__(self, 'Exception', parent.window, 0, (Gtk.DialogFlags.DESTROY_WITH_PARENT), (Gtk.STOCK_OK, Gtk.ResponseType.OK))
        self.set_border_width(10)
        box = self.get_content_area()
        label = Gtk.Label('Starterator has reached an error.')
        box.add(label)
        self.show_all()

class StarteratorFinishedDialog(Gtk.Dialog):

    def __init__(self, parent, report_file):
        Gtk.Dialog.__init__(self, 'Starterator is Done', parent, 0)
        self.set_border_width(10)
        box = self.get_content_area()
        label = Gtk.Label('Starterator has completed.')
        box.pack_start(label, False, False, 0)
        link = Gtk.Button('Link to Report')
        link.connect('clicked', self.on_link_button_clicked, report_file)
        box.pack_start(link, False, False, 0)
        again_label = Gtk.Label('Do you want to run another report?')
        box.pack_start(again_label, False, False, 0)
        self.add_button('Yes', Gtk.ResponseType.YES)
        self.add_button('No, quit Starterator', Gtk.ResponseType.NO) 
        self.show_all()

    def on_link_button_clicked(self, button, link):
        os.system('xdg-open ' + "\""+link+"\"")

class DatabaseInfoDialog(Gtk.Dialog):
    def __init__(self, parent, db_server, db_name, db_user):
        Gtk.Dialog.__init__(self, 'Database Information', parent.window, 0, (Gtk.STOCK_OK, Gtk.ResponseType.OK))
        self.set_border_width(10)
        self.info = {'db_server' : db_server,
                     'db_name' : db_name,
                     'db_user': db_user}
        vbox = self.get_content_area()
        db_info_label = Gtk.Label('Could not connect to database, please re-enter information')
        vbox.pack_start(db_info_label, False, False, 0)

        hbox = Gtk.Box(spacing= 6)
        db_server_label = Gtk.Label("Database Server")
        db_server_entry = Gtk.Entry()
        db_server_entry.set_text(db_server)
        db_server_entry.connect("changed", self.on_entry_changed, 'db_server')
        hbox.pack_start(db_server_label, False, False, 0)
        hbox.pack_start(db_server_entry, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)

        hbox = Gtk.Box(spacing= 6)
        db_name_label = Gtk.Label("Database Name")
        db_name_entry = Gtk.Entry()
        db_name_entry.set_text(db_name)
        db_name_entry.connect("changed", self.on_entry_changed, 'db_name')
        hbox.pack_start(db_name_label, False, False, 0)
        hbox.pack_start(db_name_entry, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)

        hbox = Gtk.Box(spacing= 6)
        db_user_label = Gtk.Label("Database User")
        db_user_entry = Gtk.Entry()
        db_user_entry.set_text(db_user)
        db_user_entry.connect("changed", self.on_entry_changed, 'db_user')
        hbox.pack_start(db_user_label, False, False, 0)
        hbox.pack_start(db_user_entry, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)

        hbox = Gtk.Box(spacing= 6)
        db_pass_label = Gtk.Label("Database Password")
        db_pass_entry = Gtk.Entry()
        db_pass_entry.set_visibility(False)
        db_pass_entry.connect("changed", self.on_entry_changed, 'db_password')
        hbox.pack_start(db_pass_label, False, False, 0)
        hbox.pack_start(db_pass_entry, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)
        self.show_all()


    def on_entry_changed(self, entry, name):
        self.info[name] = entry.get_text()

class StarteratorProgressDialog(Gtk.Dialog):
    def __init__(self, parent):
        Gtk.Dialog.__init__(self, 'Starterator Progress', parent.window, 0)
        vbox = self.get_content_area()
        self.progress_label = Gtk.Label('Starterator is starting')
        self.progress_bar = Gtk.ProgressBar()
        vbox.pack_start(self.progress_label, True, True, 0)
        vbox.pack_start(self.progress_bar, True, True, 0)
        self.show_all()

    def update_starterator(self, label_text, progress_amount):
        self.progress_label.set_text(label_text)
        self.progress_bar.set_fraction(progress_amount)

class StarteratorThread(threading.Thread):
    def __init__(self, parent, db, config, info):
        threading.Thread.__init__(self)
        # self.stop = False
        self.parent = parent
        self.db = db
        self.config_info = config
        self.info = info
        self.final_file = None
        self.setDaemon(True)
        self.stop_thread = threading.Event()

    def run(self):
        self.starterate()

    def starterate(self):
        try:
            self.final_file = starterates.starterate(self.db, self.config_info, self.info, 
                gui=self, event=self.stop_thread)
            # print self.final_file
        except:
            if self.stop_thread.is_set():
                return
            self.db.close()
            self.stop = True
            Gdk.threads_enter()
            dialog = Gtk.MessageDialog(self.parent, 0, Gtk.MessageType.ERROR,
                Gtk.ButtonsType.CANCEL, "Starterator has encountered an error")
            dialog.format_secondary_text(
            "Please try again.")
            dialog.run()
            dialog.destroy()
            Gdk.threads_leave()
            # self.stop = True
            raise
 
        else:
            self.db.close()
            if self.stop_thread.is_set():
                # clean up files?
                # show dialog
                return
            Gdk.threads_enter()
            done_dialog = StarteratorFinishedDialog(self.parent, self.final_file)
            response = done_dialog.run()
            if response == Gtk.ResponseType.NO:
                Gtk.main_quit()
            done_dialog.destroy()
            Gdk.threads_leave()
            # self.stop = True

    def cancel(self):
        self.stop_thread.set()

    def update(self, text, amount):
        print('update', text, amount)
        Gdk.threads_enter()
        self.parent.update_starterator(text, amount)
        Gdk.threads_leave()

class PreferencesDialog(Gtk.Dialog):
    def __init__(self, parent, config):
        Gtk.Dialog.__init__(self, "Preferences", parent.window, 0,  
             (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
             Gtk.STOCK_OK, Gtk.ResponseType.APPLY))
        self.set_border_width(10)
        self.config_info = config
        box = self.get_content_area()
        self.box = box
        gtknotebook = Gtk.Notebook()
        self.box.add(gtknotebook)
        database_box = self.database_tab()
        gtknotebook.append_page(database_box, Gtk.Label("Database"))
        file_box = self.file_tab()
        gtknotebook.append_page(file_box, Gtk.Label("File"))
        # other_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        # gtknotebook.append_page(other_box, Gtk.Label("Other"))
        self.show_all()
  

    def database_tab(self):
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        hbox = Gtk.Box(spacing=6)
        db_server_label = Gtk.Label("Database Server")
        db_server_entry = Gtk.Entry()
        db_server_entry.set_text(self.config_info["database_server"])
        db_server_entry.connect("changed", self.on_entry_changed, 'database_server')
        hbox.pack_start(db_server_label, False, False, 0)
        hbox.pack_start(db_server_entry, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)

        hbox = Gtk.Box(spacing= 6)
        db_name_label = Gtk.Label("Database Name")
        db_name_entry = Gtk.Entry()
        db_name_entry.set_text(self.config_info["database_name"])
        db_name_entry.connect("changed", self.on_entry_changed, 'database_name')
        hbox.pack_start(db_name_label, False, False, 0)
        hbox.pack_start(db_name_entry, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)

        hbox = Gtk.Box(spacing= 6)
        db_user_label = Gtk.Label("Database User")
        db_user_entry = Gtk.Entry()
        db_user_entry.set_text(self.config_info["database_user"])
        db_user_entry.connect("changed", self.on_entry_changed, 'database_user')
        hbox.pack_start(db_user_label, False, False, 0)
        hbox.pack_start(db_user_entry, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)

        hbox = Gtk.Box(spacing= 6)
        db_pass_label = Gtk.Label("Database Password")
        db_pass_entry = Gtk.Entry()
        db_pass_entry.set_visibility(False)
        db_pass_entry.set_text(self.config_info["database_password"])
        db_pass_entry.connect("changed", self.on_entry_changed, 'database_password')
        hbox.pack_start(db_pass_label, False, False, 0)
        hbox.pack_start(db_pass_entry, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)
        return vbox

    def on_entry_changed(self, entry, name):
        self.config_info[name] = entry.get_text()
        utils.write_to_config_file(self.config_info)
    
    def file_tab(self):
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        hbox = Gtk.Box(spacing= 6)
        int_button = Gtk.Button('Choose Intermediate File Folder')
        int_entry = Gtk.Entry()
        int_entry.set_text(self.config_info["intermediate_file_dir"])
        int_entry.connect('changed', self.on_entry_changed, 'intermediate_file_dir')
        int_button.connect('clicked', self.on_folder_clicked, ['intermediate_file_dir', int_entry, "Intermediate Files"])
        hbox.pack_start(int_entry, False, False, 0 )
        hbox.pack_start(int_button, False, False, 0 )
        box.pack_start(hbox, False, False, 0)

        hbox = Gtk.Box(spacing= 6)
        final_button = Gtk.Button('Choose Report Files Folder')
        final_entry = Gtk.Entry()
        final_entry.set_text(self.config_info["final_file_dir"])
        final_entry.connect('changed', self.on_entry_changed, 'final_file_dir')
        final_button.connect('clicked', self.on_folder_clicked, ['final_file_dir', final_entry, "Report Files"])
        hbox.pack_start(final_entry, False, False, 0 )
        hbox.pack_start(final_button, False, False, 0 )
        box.pack_start(hbox, False, False, 0)

        # hbox = Gtk.Box(spacing= 6)
        # protein_button = Gtk.Button('Choose a Protein Database Folder')
        # protein_entry = Gtk.Entry()
        # protein_entry.set_text(self.config_info["protein_db"])
        # protein_entry.connect('changed', self.on_entry_changed, "protein_db")
        # protein_button.connect('clicked', self.on_folder_clicked, ['protein_db', protein_entry, "Protein Databse"])
        # hbox.pack_start(protein_entry, False, False, 0 )
        # hbox.pack_start(protein_button, False, False, 0 )
        # box.pack_start(hbox, False, False, 0)
        return box

    def on_folder_clicked(self, button, data):
        name = data[0]
        entry = data[1]
        read = data[2]
        dialog = Gtk.FileChooserDialog("Please choose a folder for %s" % read, self,
            Gtk.FileChooserAction.CREATE_FOLDER,
            (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
             Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            print("File selected: " + dialog.get_filename())
            self.config_info[name] = dialog.get_filename()
            entry.set_text(self.config_info[name])
            print(name, self.config_info[name])
            utils.write_to_config_file(self.config_info)

        dialog.destroy()

    def other_tab(self):
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        return box
