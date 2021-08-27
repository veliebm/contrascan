% Delete the lock file associated with a script.
function delete_lock_file(script_name)
    % Delete the lock file associated with a script.
    lock_name = append(script_name, '.m_LOCKFILE')
    if exist(lock_name, 'file')==2
        delete(lock_name);
    end
end
